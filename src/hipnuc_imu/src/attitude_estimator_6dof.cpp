/*======================================================================
 *  AHRS 6轴算法 (无磁力计)
 *
 *  导航系:  东北天 (ENU)      X=East, Y=North, Z=Up
 *  机体系:  右前上 (RFU)      X=Right, Y=Forward, Z=Up
 *  旋转:    ZXY 内旋          yaw(绕Z) → pitch(绕X) → roll(绕Y)
 *  四元数:  q_0=w, q_1=x, q_2=y, q_3=z
 *  欧拉角:  roll=绕y, pitch=绕X, yaw=绕Z
 *
 *  注意: 6轴无磁力计时, yaw会随时间漂移(加速度计只能修正roll和pitch)
 *======================================================================*/
#include <memory>
#include <cmath>
#include <chrono>

#include "rclcpp/rclcpp.hpp"
#include "sensor_msgs/msg/imu.hpp"
#include "geometry_msgs/msg/vector3_stamped.hpp"
#include "tf2_ros/transform_broadcaster.h"
#include "tf2_geometry_msgs/tf2_geometry_msgs.hpp"

using namespace std::chrono_literals;
// ========== 常量定义 ==========
#define DEG_TO_RAD  0.017453292519943295769236907684886f
#define RAD_TO_DEG  57.295779513082320876798154814105f
namespace imu_6dof
{

// ========== 一阶低通滤波器 (EMA) ==========
class LowPassFilter
{
public:
    LowPassFilter() : initialized_(false), value_(0.0) {}

    // cutoff_freq: 截止频率 (Hz), dt: 采样周期 (s)
    void set_params(double cutoff_freq, double dt)
    {
        if (cutoff_freq <= 0.0) {
            alpha_ = 1.0;  // 不滤波
        } else {
            double rc = 1.0 / (2.0 * M_PI * cutoff_freq);
            alpha_ = dt / (dt + rc);
            // 限制 alpha 范围，避免极端值
            alpha_ = std::max(0.01, std::min(1.0, alpha_));
        }
    }

    double update(double new_value)
    {
        if (!initialized_) {
            value_ = new_value;
            initialized_ = true;
            return value_;
        }
        value_ = alpha_ * new_value + (1.0 - alpha_) * value_;
        return value_;
    }

    double get() const { return value_; }
    void reset(double v = 0.0) { value_ = v; initialized_ = false; }

private:
    bool initialized_;
    double value_;
    double alpha_ = 1.0;
};


// ========== 6-DOF 姿态解算节点 ==========
class MahonyMadgwick6DOF : public rclcpp::Node
{
public:
    MahonyMadgwick6DOF() : Node("imu_6dof_node")
    // ROS的构造函数
    {
        // ---- 声明参数 ----
        this->declare_parameter<std::string>("imu_topic", "/IMU_data");
        this->declare_parameter<std::string>("output_topic", "/attitude");
        this->declare_parameter<int>("algorithm", 1);            // 0: Mahony, 1: Madgwick
        this->declare_parameter<double>("kp", 0.2);
        this->declare_parameter<double>("ki", 0.001);
        this->declare_parameter<double>("beta", 0.05);
        this->declare_parameter<double>("init_static_time", 2.0);

        // 低通滤波参数
        this->declare_parameter<double>("accel_lpf_freq", 1.0);  // 加速度计截止频率 Hz
        this->declare_parameter<double>("gyro_lpf_freq", 1.0);  // 陀螺仪截止频率 Hz
        this->declare_parameter<bool>("enable_lpf", true);       // 是否启用滤波

        // ---- 获取参数 ----
        this->get_parameter("imu_topic", imu_topic_);
        this->get_parameter("output_topic", output_topic_);
        this->get_parameter("algorithm", algorithm_);
        this->get_parameter("kp", kp_);
        this->get_parameter("ki", ki_);
        this->get_parameter("beta", beta_);
        this->get_parameter("init_static_time", init_static_time_);
        this->get_parameter("accel_lpf_freq", accel_lpf_freq_);
        this->get_parameter("gyro_lpf_freq", gyro_lpf_freq_);
        this->get_parameter("enable_lpf", enable_lpf_);

        // ---- ROS 接口 ----
        sub_imu_ = this->create_subscription<sensor_msgs::msg::Imu>(
            imu_topic_, rclcpp::SensorDataQoS(),
            std::bind(&MahonyMadgwick6DOF::imu_callback, this, std::placeholders::_1));

        pub_attitude_ = this->create_publisher<geometry_msgs::msg::Vector3Stamped>(output_topic_, 10);

        pub_filtered_imu_ = this->create_publisher<sensor_msgs::msg::Imu>("/imu_filtered", 10);

        tf_broadcaster_ = std::make_shared<tf2_ros::TransformBroadcaster>(this);

        // ---- 状态初始化 ----
        q0_ = 1.0; q1_ = 0.0; q2_ = 0.0; q3_ = 0.0;
        exInt_ = 0.0; eyInt_ = 0.0; ezInt_ = 0.0;
        initialized_ = false;
        init_sample_count_ = 0;
        init_ax_sum_ = init_ay_sum_ = init_az_sum_ = 0.0;
        lpf_configured_ = false;

        last_time_ = this->now();

        const char* algo_name = (algorithm_ == 0) ? "Mahony" : "Madgwick";
        RCLCPP_INFO(this->get_logger(),
            "6-DOF %s filter started | LPF %s (accel=%.1fHz, gyro=%.1fHz)",
            algo_name, enable_lpf_ ? "ON" : "OFF", accel_lpf_freq_, gyro_lpf_freq_);
    }

private:
    //接收到IMU数据后就开始执行的回调函数
    void imu_callback(const sensor_msgs::msg::Imu::SharedPtr msg)
    {
        auto now = this->now();
        double dt = (now - last_time_).seconds();
        if (dt <= 0.0 || dt > 0.1) {
            last_time_ = now;
            return;
        }

        // 首次收到数据时配置滤波器的 alpha（需要真实 dt）
        if (!lpf_configured_ && enable_lpf_) {
            lpf_ax_.set_params(accel_lpf_freq_, dt);
            lpf_ay_.set_params(accel_lpf_freq_, dt);
            lpf_az_.set_params(accel_lpf_freq_, dt);
            lpf_gx_.set_params(gyro_lpf_freq_, dt);
            lpf_gy_.set_params(gyro_lpf_freq_, dt);
            lpf_gz_.set_params(gyro_lpf_freq_, dt);
            lpf_configured_ = true;
        }

        // 读取原始数据
        double gx = msg->angular_velocity.x;
        double gy = msg->angular_velocity.y;
        double gz = msg->angular_velocity.z;
        double ax = msg->linear_acceleration.x;
        double ay = msg->linear_acceleration.y;
        double az = msg->linear_acceleration.z;

        // ---- 低通滤波 ----
        if (enable_lpf_) {
            ax = lpf_ax_.update(ax);
            ay = lpf_ay_.update(ay);
            az = lpf_az_.update(az);
            gx = lpf_gx_.update(gx);
            gy = lpf_gy_.update(gy);
            gz = lpf_gz_.update(gz);
        }

        // 发布滤波后的 IMU 数据（方便调试对比）
        auto filtered_msg = std::make_unique<sensor_msgs::msg::Imu>();
        filtered_msg->header = msg->header;
        filtered_msg->angular_velocity.x = gx;
        filtered_msg->angular_velocity.y = gy;
        filtered_msg->angular_velocity.z = gz;
        filtered_msg->linear_acceleration.x = ax;
        filtered_msg->linear_acceleration.y = ay;
        filtered_msg->linear_acceleration.z = az;
        filtered_msg->orientation_covariance = msg->orientation_covariance;
        filtered_msg->angular_velocity_covariance = msg->angular_velocity_covariance;
        filtered_msg->linear_acceleration_covariance = msg->linear_acceleration_covariance;
        pub_filtered_imu_->publish(std::move(filtered_msg));

        // ---- 初始对准 ----
        if (!initialized_) {
            init_ax_sum_ += ax;
            init_ay_sum_ += ay;
            init_az_sum_ += az;
            init_sample_count_++;
            double elapsed = init_sample_count_ * dt;
            if (elapsed >= init_static_time_) {
                double avg_ax = init_ax_sum_ / init_sample_count_;
                double avg_ay = init_ay_sum_ / init_sample_count_;
                double avg_az = init_az_sum_ / init_sample_count_;
                init_attitude(avg_ax, avg_ay, avg_az);
                initialized_ = true;
                RCLCPP_INFO(this->get_logger(), "6-DOF initialized (%d samples)", init_sample_count_);
                
            }
            last_time_ = now;
            return;
        }

        // ---- 滤波更新 ----
        if (algorithm_ == 0) {
            mahony_update(gx, gy, gz, ax, ay, az, dt);
        } else {
            madgwick_update(gx, gy, gz, ax, ay, az, dt);
        }

        publish_attitude(now);
        broadcast_tf(now);
        last_time_ = now;
    }

    // 初始化 由静止加速度计计算初始 Roll 和 Pitch
    void init_attitude(double ax, double ay, double az)
    {
        double norm = sqrt(ax * ax + ay * ay + az * az);
        if (norm < 1e-6) return;
        ax /= norm; ay /= norm; az /= norm;

        double roll  = atan2(-ax, az);
        double pitch = atan2(ay, sqrt(ax * ax + az * az));
        double yaw = 0.0;
        
        RCLCPP_INFO(this->get_logger(), 
        "Initial attitude: roll=%.2f°, pitch=%.2f°, yaw=%.2f°",
        roll * RAD_TO_DEG, pitch * RAD_TO_DEG, yaw * RAD_TO_DEG);

        // Z-X-Y 欧拉角转四元数
        double cy = cos(yaw * 0.5), sy = sin(yaw * 0.5);
        double cp = cos(pitch * 0.5), sp = sin(pitch * 0.5);
        double cr = cos(roll * 0.5), sr = sin(roll * 0.5);
        q0_ = cy * cp * cr - sy * sp * sr;
        q1_ = cy * sp * cr - sy * cp * sr;
        q2_ = cy * cp * sr + sy * sp * cr;
        q3_ = sy * cp * cr + cy * sp * sr;

        exInt_ = eyInt_ = ezInt_ = 0.0;
    }

    // Mahony 滤波 PIPI控制
    void mahony_update(double gx, double gy, double gz,
                       double ax, double ay, double az, double dt)
    {
        double norm = sqrt(ax * ax + ay * ay + az * az);
        if (norm < 1e-6) return;
        ax /= norm; ay /= norm; az /= norm;

        double vx = 2.0 * (q1_ * q3_ - q0_ * q2_);
        double vy = 2.0 * (q0_ * q1_ + q2_ * q3_);
        double vz = q0_ * q0_ - q1_ * q1_ - q2_ * q2_ + q3_ * q3_;

        double ex = ay * vz - az * vy;
        double ey = az * vx - ax * vz;
        double ez = ax * vy - ay * vx;

        exInt_ += ex * ki_ * dt;
        eyInt_ += ey * ki_ * dt;
        ezInt_ += ez * ki_ * dt;

        const double INTEGRAL_LIMIT = 0.3;
        exInt_ = std::max(-INTEGRAL_LIMIT, std::min(INTEGRAL_LIMIT, exInt_));
        eyInt_ = std::max(-INTEGRAL_LIMIT, std::min(INTEGRAL_LIMIT, eyInt_));
        ezInt_ = std::max(-INTEGRAL_LIMIT, std::min(INTEGRAL_LIMIT, ezInt_));

        gx += kp_ * ex + exInt_;
        gy += kp_ * ey + eyInt_;
        gz += kp_ * ez + ezInt_;

        double q0_dot = 0.5 * (-q1_ * gx - q2_ * gy - q3_ * gz);
        double q1_dot = 0.5 * ( q0_ * gx + q2_ * gz - q3_ * gy);
        double q2_dot = 0.5 * ( q0_ * gy - q1_ * gz + q3_ * gx);
        double q3_dot = 0.5 * ( q0_ * gz + q1_ * gy - q2_ * gx);
        q0_ += q0_dot * dt;
        q1_ += q1_dot * dt;
        q2_ += q2_dot * dt;
        q3_ += q3_dot * dt;
        normalize_quat();
    }

    // Madgwick 滤波
    void madgwick_update(double gx, double gy, double gz,
                         double ax, double ay, double az, double dt)
    {
        double norm = sqrt(ax * ax + ay * ay + az * az);
        if (norm < 1e-6) return;
        ax /= norm; ay /= norm; az /= norm;

        double q0 = q0_, q1 = q1_, q2 = q2_, q3 = q3_;

        double f1 = 2.0 * (q1 * q3 - q0 * q2) - ax;
        double f2 = 2.0 * (q0 * q1 + q2 * q3) - ay;
        double f3 = 2.0 * (0.5 - q1 * q1 - q2 * q2) - az;

        double grad0 = -2 * q2 * f1 + 2 * q1 * f2;
        double grad1 =  2 * q3 * f1 + 2 * q0 * f2 - 4 * q1 * f3;
        double grad2 = -2 * q0 * f1 + 2 * q3 * f2 - 4 * q2 * f3;
        double grad3 =  2 * q1 * f1 + 2 * q2 * f2;

        norm = sqrt(grad0 * grad0 + grad1 * grad1 + grad2 * grad2 + grad3 * grad3);
        if (norm > 1e-6) {
            grad0 /= norm; grad1 /= norm; grad2 /= norm; grad3 /= norm;
        }

        double qDot0 = 0.5 * (-q1 * gx - q2 * gy - q3 * gz);
        double qDot1 = 0.5 * ( q0 * gx + q2 * gz - q3 * gy);
        double qDot2 = 0.5 * ( q0 * gy - q1 * gz + q3 * gx);
        double qDot3 = 0.5 * ( q0 * gz + q1 * gy - q2 * gx);

        qDot0 -= beta_ * grad0;
        qDot1 -= beta_ * grad1;
        qDot2 -= beta_ * grad2;
        qDot3 -= beta_ * grad3;

        q0_ += qDot0 * dt;
        q1_ += qDot1 * dt;
        q2_ += qDot2 * dt;
        q3_ += qDot3 * dt;
        normalize_quat();
    }

    void normalize_quat()
    {
        double norm = sqrt(q0_ * q0_ + q1_ * q1_ + q2_ * q2_ + q3_ * q3_);
        if (norm > 1e-6) {
            q0_ /= norm; q1_ /= norm; q2_ /= norm; q3_ /= norm;
        }
    }

    void quat_to_euler(double w, double x, double y, double z,
                       double& roll, double& pitch, double& yaw)
    {
        pitch = asin(2.0 * (w * x + y * z));
        roll  = -atan2(2.0 * (x * z - w * y), w * w - x * x - y * y + z * z);
        yaw   = -atan2(2.0 * (x * y - w * z), w * w - x * x + y * y - z * z);
    }
    // ---- ROS ----
    void publish_attitude(rclcpp::Time stamp)
    {
        geometry_msgs::msg::Vector3Stamped euler_msg;
        euler_msg.header.stamp = stamp;
        euler_msg.header.frame_id = "imu_link";

        double roll, pitch, yaw;

        //这一步也可以直接输出四元数了
        quat_to_euler(q0_, q1_, q2_, q3_, roll, pitch, yaw);
        euler_msg.vector.x = roll;
        euler_msg.vector.y = pitch;
        euler_msg.vector.z = yaw;
        pub_attitude_->publish(euler_msg);

        static int cnt = 0;
        if (++cnt % 100 == 0) {
            RCLCPP_INFO(this->get_logger(),
                        "Roll: %6.2f deg, Pitch: %6.2f deg, Yaw: %6.2f deg",
                        roll * 180.0 / M_PI, pitch * 180.0 / M_PI, yaw * 180.0 / M_PI);
        }
    }
    // ---- ROS ----
    void broadcast_tf(rclcpp::Time stamp)
    {
        geometry_msgs::msg::TransformStamped tf;
        tf.header.stamp = stamp;
        tf.header.frame_id = "world";
        tf.child_frame_id = "imu_attitude";
        tf.transform.rotation.w = q0_;
        tf.transform.rotation.x = q1_;
        tf.transform.rotation.y = q2_;
        tf.transform.rotation.z = q3_;
        tf_broadcaster_->sendTransform(tf);
    }

    // ---- 参数 ----
    std::string imu_topic_, output_topic_;
    int algorithm_;
    double kp_, ki_, beta_;
    double init_static_time_;
    double accel_lpf_freq_, gyro_lpf_freq_;
    bool enable_lpf_;

    // ---- 低通滤波器 ----
    bool lpf_configured_;
    LowPassFilter lpf_ax_, lpf_ay_, lpf_az_;
    LowPassFilter lpf_gx_, lpf_gy_, lpf_gz_;

    // ---- 状态 ----
    double q0_, q1_, q2_, q3_;
    double exInt_, eyInt_, ezInt_;

    // ---- 初始化 ----
    bool initialized_;
    int init_sample_count_;
    double init_ax_sum_, init_ay_sum_, init_az_sum_;

    // ---- ROS ----
    rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr sub_imu_;
    rclcpp::Publisher<geometry_msgs::msg::Vector3Stamped>::SharedPtr pub_attitude_;
    rclcpp::Publisher<sensor_msgs::msg::Imu>::SharedPtr pub_filtered_imu_;
    std::shared_ptr<tf2_ros::TransformBroadcaster> tf_broadcaster_;
    rclcpp::Time last_time_;
};

} // namespace imu_6dof

int main(int argc, char** argv)
{
    rclcpp::init(argc, argv);
    auto node = std::make_shared<imu_6dof::MahonyMadgwick6DOF>();
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}
