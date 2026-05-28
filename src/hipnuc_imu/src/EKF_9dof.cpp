//yaw角没有实现
/*======================================================================
 *  EKF 9轴 AHRS (含磁力计)
 *
 *  导航系:  东北天 (ENU)      X=East, Y=North, Z=Up
 *  机体系:  右前上 (RFU)      X=Right, Y=Forward, Z=Up
 *  旋转:    ZXY 内旋          yaw(绕Z) → pitch(绕X) → roll(绕Y)
 *  四元数:  q_0=w, q_1=x, q_2=y, q_3=z (Hamilton)
 *  欧拉角:  roll=绕Y, pitch=绕X, yaw=绕Z
 *
 *  状态向量 (7D): [qw, qx, qy, qz, bias_x, bias_y, bias_z]
 *  预测量: 陀螺仪驱动四元数运动学
 *  观测量: 加速度计(重力方向) + 磁力计(地磁方向)
 *
 *  欧拉角转四元数(ZXY): q = qz(yaw) ⊗ qx(pitch) ⊗ qy(roll)
 *    q0 = cy*cp*cr - sy*sp*sr
 *    q1 = cy*sp*cr - sy*cp*sr
 *    q2 = cy*cp*sr + sy*sp*cr
 *    q3 = sy*cp*cr + cy*sp*sr
 *
 *  四元数转欧拉角(ZXY):
 *    pitch = arcsin(2*(w*x + y*z))
 *    roll  = -atan2(2*(x*z - w*y), w^2 - x^2 - y^2 + z^2)
 *    yaw   = -atan2(2*(x*y - w*z), w^2 - x^2 + y^2 - z^2)
 * // 方向余弦矩阵 C_b2n (从机体坐标系到导航坐标系)
// 使用四元数 q = [q0, q1, q2, q3] 表示
//
// |  q0²+q1²-q2²-q3²    2(q1q2-q0q3)      2(q1q3+q0q2)   |
// |                                                       |
// |  2(q1q2+q0q3)       q0²-q1²+q2²-q3²   2(q2q3-q0q1)   |
// |                                                       |
// |  2(q1q3-q0q2)       2(q2q3+q0q1)      q0²-q1²-q2²+q3²|
 *======================================================================*/
#include <memory>
#include <cmath>
#include <chrono>
#include <cstdint>

#include "rclcpp/rclcpp.hpp"
#include "sensor_msgs/msg/imu.hpp"
#include "sensor_msgs/msg/magnetic_field.hpp"
#include "geometry_msgs/msg/vector3_stamped.hpp"
#include "tf2_ros/transform_broadcaster.h"
#include "tf2_geometry_msgs/tf2_geometry_msgs.hpp"

using namespace std::chrono_literals;



#define DEG_TO_RAD  0.017453292519943295769236907684886
#define RAD_TO_DEG  57.295779513082320876798154814105

namespace imu_9dof
{

// ========== 快速开方倒数 ==========
static double inv_sqrt(double x)
{
    double half = 0.5 * x;
    int64_t i  = *(int64_t*)&x;
    i = 0x5FE6EB50C7B537A9 - (i >> 1);
    x = *(double*)&i;
    x = x * (1.5 - half * x * x);
    x = x * (1.5 - half * x * x);
    return x;
}

// ========== 一阶低通滤波器 ==========
class LowPassFilter
{
public:
    LowPassFilter() : initialized_(false), value_(0.0) {}

    void set_params(double cutoff_freq, double dt)
    {
        if (cutoff_freq <= 0.0) {
            alpha_ = 1.0;
        } else {
            double rc = 1.0 / (2.0 * M_PI * cutoff_freq);
            alpha_ = dt / (dt + rc);
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

// ========== 矩阵运算辅助 ==========
static void mat_mul_7x7(const double A[7][7], const double B[7][7], double C[7][7])
{
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 7; j++) {
            double sum = 0.0;
            for (int k = 0; k < 7; k++) sum += A[i][k] * B[k][j];
            C[i][j] = sum;
        }
}

static void mat_mul_7x7_AT(const double A[7][7], const double B[7][7], double C[7][7])
{
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 7; j++) {
            double sum = 0.0;
            for (int k = 0; k < 7; k++) sum += A[i][k] * B[j][k];
            C[i][j] = sum;
        }
}

static void mat_mul_3x7_7x7(const double A[3][7], const double B[7][7], double C[3][7])
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 7; j++) {
            double sum = 0.0;
            for (int k = 0; k < 7; k++) sum += A[i][k] * B[k][j];
            C[i][j] = sum;
        }
}

static bool mat_inv_3x3(const double S[3][3], double invS[3][3])
{
    double det = S[0][0] * (S[1][1] * S[2][2] - S[1][2] * S[2][1])
               - S[0][1] * (S[1][0] * S[2][2] - S[1][2] * S[2][0])
               + S[0][2] * (S[1][0] * S[2][1] - S[1][1] * S[2][0]);
    if (std::fabs(det) < 1e-15) return false;

    double invDet = 1.0 / det;
    invS[0][0] = (S[1][1] * S[2][2] - S[1][2] * S[2][1]) * invDet;
    invS[0][1] = (S[0][2] * S[2][1] - S[0][1] * S[2][2]) * invDet;
    invS[0][2] = (S[0][1] * S[1][2] - S[0][2] * S[1][1]) * invDet;
    invS[1][0] = (S[1][2] * S[2][0] - S[1][0] * S[2][2]) * invDet;
    invS[1][1] = (S[0][0] * S[2][2] - S[0][2] * S[2][0]) * invDet;
    invS[1][2] = (S[0][2] * S[1][0] - S[0][0] * S[1][2]) * invDet;
    invS[2][0] = (S[1][0] * S[2][1] - S[1][1] * S[2][0]) * invDet;
    invS[2][1] = (S[0][1] * S[2][0] - S[0][0] * S[2][1]) * invDet;
    invS[2][2] = (S[0][0] * S[1][1] - S[0][1] * S[1][0]) * invDet;
    return true;
}


// ========== 9-DOF EKF 姿态解算节点 ==========
class EKF9DOF : public rclcpp::Node
{
public:
    EKF9DOF() : Node("ekf_9dof_node")
    {
        // ---- 参数声明 ----
        this->declare_parameter<std::string>("imu_topic", "/IMU_data");
        this->declare_parameter<std::string>("mag_topic", "/magnetic_data");
        this->declare_parameter<std::string>("euler_topic", "/euler_data");
        this->declare_parameter<std::string>("output_topic", "/attitude");
        this->declare_parameter<double>("init_static_time", 2.0);

        // EKF噪声参数
        this->declare_parameter<double>("sigma_g", 0.02);       // 陀螺仪角速度噪声 rad/s/√Hz
        this->declare_parameter<double>("sigma_b", 1.0e-4);     // 零偏随机游走 rad/s/√s/√Hz
        this->declare_parameter<double>("sigma_a", 0.01);       // 加速度计量测噪声 (归一化)
        this->declare_parameter<double>("sigma_m", 0.05);       // 磁力计量测噪声 (归一化)

        // 低通滤波参数
        this->declare_parameter<double>("accel_lpf_freq", 30.0);
        this->declare_parameter<double>("gyro_lpf_freq", 30.0);
        this->declare_parameter<double>("mag_lpf_freq", 10.0);   // 磁力计截止频率 (通常数据率更低)
        this->declare_parameter<bool>("enable_lpf", false);

        // 量测门控
        this->declare_parameter<double>("expected_g_norm", 9.81); // 期望重力幅值
        this->declare_parameter<double>("accel_gate", 2.0);       // 外部加速度阈值
        this->declare_parameter<double>("mag_gate", 0.3);         // 磁干扰阈值 (归一化幅值偏差)

        // ---- 获取参数 ----
        this->get_parameter("imu_topic", imu_topic_);
        this->get_parameter("mag_topic", mag_topic_);
        this->get_parameter("euler_topic", euler_topic_);
        this->get_parameter("output_topic", output_topic_);
        this->get_parameter("init_static_time", init_static_time_);
        this->get_parameter("sigma_g", sigma_g_);
        this->get_parameter("sigma_b", sigma_b_);
        this->get_parameter("sigma_a", sigma_a_);
        this->get_parameter("sigma_m", sigma_m_);
        this->get_parameter("accel_lpf_freq", accel_lpf_freq_);
        this->get_parameter("gyro_lpf_freq", gyro_lpf_freq_);
        this->get_parameter("mag_lpf_freq", mag_lpf_freq_);
        this->get_parameter("enable_lpf", enable_lpf_);
        this->get_parameter("expected_g_norm", expected_g_norm_);
        this->get_parameter("accel_gate", accel_gate_);
        this->get_parameter("mag_gate", mag_gate_);

        // ---- ROS 接口 ----
        // 使用回调组实现多话题同步，或简单缓存磁力计数据
        sub_imu_ = this->create_subscription<sensor_msgs::msg::Imu>(
            imu_topic_, rclcpp::SensorDataQoS(),
            std::bind(&EKF9DOF::imu_callback, this, std::placeholders::_1));

        sub_mag_ = this->create_subscription<sensor_msgs::msg::MagneticField>(
            mag_topic_, rclcpp::SensorDataQoS(),
            std::bind(&EKF9DOF::mag_callback, this, std::placeholders::_1));

        sub_euler_ = this->create_subscription<geometry_msgs::msg::Vector3Stamped>(
            euler_topic_, rclcpp::SensorDataQoS(),
            std::bind(&EKF9DOF::euler_callback, this, std::placeholders::_1));

        pub_attitude_ = this->create_publisher<geometry_msgs::msg::Vector3Stamped>(output_topic_, 10);
        pub_filtered_imu_ = this->create_publisher<sensor_msgs::msg::Imu>("/imu_filtered", 10);
        tf_broadcaster_ = std::make_shared<tf2_ros::TransformBroadcaster>(this);

        // ---- 状态初始化 ----
        x_[0] = 1.0; x_[1] = 0.0; x_[2] = 0.0; x_[3] = 0.0;
        x_[4] = 0.0; x_[5] = 0.0; x_[6] = 0.0;
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < 7; j++)
                P_[i][j] = (i == j) ? ((i < 4) ? 0.001 : 0.0001) : 0.0;

        ekf_inited_ = false;
        initialized_ = false;
        mag_received_ = false;
        init_sample_count_ = 0;
        init_ax_sum_ = init_ay_sum_ = init_az_sum_ = 0.0;
        init_mx_sum_ = init_my_sum_ = init_mz_sum_ = 0.0;
        lpf_configured_ = false;

        // 磁力计参考向量 (在导航系ENU下，初始化时确定)
        mag_ref_[0] = 0.0;  // East=0
        mag_ref_[1] = 0.0;  // North
        mag_ref_[2] = 0.0;  // Up (inclination补偿后)

        last_time_ = this->now();

        RCLCPP_INFO(this->get_logger(),
            "9-DOF EKF started | sigma_g=%.3f sigma_b=%.1e sigma_a=%.3f sigma_m=%.3f | "
            "g_norm=%.2f accel_gate=%.2f mag_gate=%.2f | LPF %s",
            sigma_g_, sigma_b_, sigma_a_, sigma_m_,
            expected_g_norm_, accel_gate_, mag_gate_,
            enable_lpf_ ? "ON" : "OFF");
    }

private:
    // ===================== 磁力计回调 (缓存最新数据) =====================
    void mag_callback(const sensor_msgs::msg::MagneticField::SharedPtr msg)
    {
        raw_mx_ = msg->magnetic_field.x;
        raw_my_ = msg->magnetic_field.y;
        raw_mz_ = msg->magnetic_field.z;
        mag_received_ = true;
    }

    // 缓存硬件欧拉角 (用于设置初始yaw)
    void euler_callback(const geometry_msgs::msg::Vector3Stamped::SharedPtr msg)
    {
        hw_roll_  = msg->vector.x;
        hw_pitch_ = msg->vector.y;
        hw_yaw_   = msg->vector.z;
        euler_received_ = true;
    }

    // ===================== IMU 回调 (主滤波循环) =====================
    void imu_callback(const sensor_msgs::msg::Imu::SharedPtr msg)
    {
        auto now = this->now();
        double dt = (now - last_time_).seconds();
        if (dt <= 0.0 || dt > 0.1) {
            last_time_ = now;
            return;
        }

        // 首次配置滤波器
        if (!lpf_configured_ && enable_lpf_) {
            lpf_ax_.set_params(accel_lpf_freq_, dt);
            lpf_ay_.set_params(accel_lpf_freq_, dt);
            lpf_az_.set_params(accel_lpf_freq_, dt);
            lpf_gx_.set_params(gyro_lpf_freq_, dt);
            lpf_gy_.set_params(gyro_lpf_freq_, dt);
            lpf_gz_.set_params(gyro_lpf_freq_, dt);
            lpf_mx_.set_params(mag_lpf_freq_, dt);
            lpf_my_.set_params(mag_lpf_freq_, dt);
            lpf_mz_.set_params(mag_lpf_freq_, dt);
            lpf_configured_ = true;
        }

        // 读取 IMU 数据
        double gx = msg->angular_velocity.x;
        double gy = msg->angular_velocity.y;
        double gz = msg->angular_velocity.z;
        double ax = msg->linear_acceleration.x;
        double ay = msg->linear_acceleration.y;
        double az = msg->linear_acceleration.z;

        // 低通滤波
        if (enable_lpf_) {
            ax = lpf_ax_.update(ax);
            ay = lpf_ay_.update(ay);
            az = lpf_az_.update(az);
            gx = lpf_gx_.update(gx);
            gy = lpf_gy_.update(gy);
            gz = lpf_gz_.update(gz);
        }

        // 发布滤波后的 IMU
        auto filtered_msg = std::make_unique<sensor_msgs::msg::Imu>();
        filtered_msg->header = msg->header;
        filtered_msg->angular_velocity.x = gx;
        filtered_msg->angular_velocity.y = gy;
        filtered_msg->angular_velocity.z = gz;
        filtered_msg->linear_acceleration.x = ax;
        filtered_msg->linear_acceleration.y = ay;
        filtered_msg->linear_acceleration.z = az;
        pub_filtered_imu_->publish(std::move(filtered_msg));

        // 磁力计滤波与缓存
        double mx = 0.0, my = 0.0, mz = 0.0;
        bool has_mag = false;
        if (mag_received_) {
            mx = enable_lpf_ ? lpf_mx_.update(raw_mx_) : raw_mx_;
            my = enable_lpf_ ? lpf_my_.update(raw_my_) : raw_my_;
            mz = enable_lpf_ ? lpf_mz_.update(raw_mz_) : raw_mz_;
            has_mag = true;
        }

        // ---- 初始对准 (需要硬件欧拉角 + 磁力计) ----
        if (!initialized_) {
            if (!euler_received_ || !has_mag) {
                last_time_ = now;
                return;  // 等待硬件欧拉角数据 + 磁力计
            }
            init_ax_sum_ += ax;  init_ay_sum_ += ay;  init_az_sum_ += az;
            init_mx_sum_ += mx;  init_my_sum_ += my;  init_mz_sum_ += mz;
            init_sample_count_++;
            if (init_sample_count_ * dt >= init_static_time_) {
                double avg_ax = init_ax_sum_ / init_sample_count_;
                double avg_ay = init_ay_sum_ / init_sample_count_;
                double avg_az = init_az_sum_ / init_sample_count_;
                double avg_mx = init_mx_sum_ / init_sample_count_;
                double avg_my = init_my_sum_ / init_sample_count_;
                double avg_mz = init_mz_sum_ / init_sample_count_;
                // 使用硬件yaw而不是磁力计推导的yaw
                init_attitude(avg_ax, avg_ay, avg_az, avg_mx, avg_my, avg_mz);
                initialized_ = true;
                RCLCPP_INFO(this->get_logger(),
                    "EKF initialized (%d samples) | hw_yaw=%.2f° | mag_ref=[%.3f %.3f %.3f]",
                    init_sample_count_, hw_yaw_ * RAD_TO_DEG,
                    mag_ref_[0], mag_ref_[1], mag_ref_[2]);
            }
            last_time_ = now;
            return;
        }

        // ---- EKF 更新 ----
        ekf_predict(gx, gy, gz, dt);

        // 顺序更新: 先加速度计，后磁力计
        ekf_update_accel(ax, ay, az);
        if (has_mag) {
            ekf_update_mag(mx, my, mz);
        }

        publish_attitude(now);
        broadcast_tf(now);
        last_time_ = now;
    }

    // ===================== 初始对准 (加速度计 + 硬件yaw + 磁力计参考) =====================
    void init_attitude(double ax, double ay, double az,
                       double mx, double my, double mz)
    {
        double norm_a = std::sqrt(ax * ax + ay * ay + az * az);
        double norm_m = std::sqrt(mx * mx + my * my + mz * mz);
        if (norm_a < 1e-6 || norm_m < 1e-6) return;

        ax /= norm_a; ay /= norm_a; az /= norm_a;
        mx /= norm_m; my /= norm_m; mz /= norm_m;

        // 由加速度计计算 roll 和 pitch
        double roll  = std::atan2(-ax, az);
        double pitch = std::atan2(ay, std::sqrt(ax * ax + az * az));
        double cp = cos(pitch), sp = sin(pitch);
        double cr = cos(roll),  sr = sin(roll);
        double Xh = cr * mx + sr * mz;
        double Yh = sp * sr * mx + cp * my - sp * cr * mz;
        double yaw = atan2(Xh, Yh);
        RCLCPP_INFO(this->get_logger(),
            "Init: roll=%.2f° pitch=%.2f° yaw=%.2f°",
            roll * RAD_TO_DEG, pitch * RAD_TO_DEG, yaw * RAD_TO_DEG);

        // ZXY 欧拉角 → 四元数: q = qz(yaw) ⊗ qx(pitch) ⊗ qy(roll)
        double cy = std::cos(yaw * 0.5),   sy = std::sin(yaw * 0.5);
        double cp_ = std::cos(pitch * 0.5), sp_ = std::sin(pitch * 0.5);
        double cr_ = std::cos(roll * 0.5),  sr_ = std::sin(roll * 0.5);

        x_[0] = cy * cp_ * cr_ - sy * sp_ * sr_;
        x_[1] = cy * sp_ * cr_ - sy * cp_ * sr_;
        x_[2] = cy * cp_ * sr_ + sy * sp_ * cr_;
        x_[3] = sy * cp_ * cr_ + cy * sp_ * sr_;
        x_[4] = 0.0; x_[5] = 0.0; x_[6] = 0.0;

        // 计算磁力计参考向量 (导航系 ENU)
        double qw = x_[0], qx = x_[1], qy = x_[2], qz = x_[3];
        double mn_x = (1.0 - 2.0*(qy*qy + qz*qz))*mx + 2.0*(qx*qy + qw*qz)*my + 2.0*(qx*qz - qw*qy)*mz;
        double mn_y = 2.0*(qx*qy - qw*qz)*mx + (1.0 - 2.0*(qx*qx + qz*qz))*my + 2.0*(qy*qz + qw*qx)*mz;
        double mn_z = 2.0*(qx*qz + qw*qy)*mx + 2.0*(qy*qz - qw*qx)*my + (1.0 - 2.0*(qx*qx + qy*qy))*mz;

        // East=0, 保留北向和天向
        double north_norm = std::sqrt(mn_x * mn_x + mn_y * mn_y);
        mag_ref_[0] = 0.0;
        mag_ref_[1] = north_norm;
        mag_ref_[2] = mn_z;

        // 初始化协方差
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < 7; j++)
                P_[i][j] = (i == j) ? ((i < 4) ? 0.001 : 0.0001) : 0.0;

        ekf_inited_ = true;
    }

    // ===================== EKF 预测 (与6轴相同) =====================
    void ekf_predict(double gx, double gy, double gz, double dt)
    {
        if (!ekf_inited_) return;

        double qw = x_[0], qx = x_[1], qy = x_[2], qz = x_[3];
        double bx = x_[4], by = x_[5], bz = x_[6];

        // 零偏补偿
        double wx = gx - bx;
        double wy = gy - by;
        double wz = gz - bz;

        // 状态预测 (Euler积分 + 归一化)
        double half_dt = 0.5 * dt;
        double qwp = qw + half_dt * (-qx * wx - qy * wy - qz * wz);
        double qxp = qx + half_dt * ( qw * wx + qy * wz - qz * wy);
        double qyp = qy + half_dt * ( qw * wy - qx * wz + qz * wx);
        double qzp = qz + half_dt * ( qw * wz + qx * wy - qy * wx);
        double qn = inv_sqrt(qwp * qwp + qxp * qxp + qyp * qyp + qzp * qzp);
        qwp *= qn; qxp *= qn; qyp *= qn; qzp *= qn;

        double bxp = bx, byp = by, bzp = bz;

        // 状态转移雅可比 F (7x7)
        double F[7][7];
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < 7; j++)
                F[i][j] = (i == j) ? 1.0 : 0.0;

        F[0][1] = -half_dt * wx;  F[0][2] = -half_dt * wy;  F[0][3] = -half_dt * wz;
        F[1][0] =  half_dt * wx;  F[1][2] =  half_dt * wz;  F[1][3] = -half_dt * wy;
        F[2][0] =  half_dt * wy;  F[2][1] = -half_dt * wz;  F[2][3] =  half_dt * wx;
        F[3][0] =  half_dt * wz;  F[3][1] =  half_dt * wy;  F[3][2] = -half_dt * wx;

        F[0][4] =  half_dt * qx;  F[0][5] =  half_dt * qy;  F[0][6] =  half_dt * qz;
        F[1][4] = -half_dt * qw;  F[1][5] =  half_dt * qz;  F[1][6] = -half_dt * qy;
        F[2][4] = -half_dt * qz;  F[2][5] = -half_dt * qw;  F[2][6] =  half_dt * qx;
        F[3][4] =  half_dt * qy;  F[3][5] = -half_dt * qx;  F[3][6] = -half_dt * qw;

        // 过程噪声协方差 Q
        double Q[7][7];
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < 7; j++)
                Q[i][j] = 0.0;

        double a = sigma_g_ * sigma_g_ * half_dt * half_dt;
        Q[0][0] = a * (1.0 - qw*qw); Q[0][1] = a * (-qw*qx); Q[0][2] = a * (-qw*qy); Q[0][3] = a * (-qw*qz);
        Q[1][0] = a * (-qw*qx);      Q[1][1] = a * (1.0 - qx*qx); Q[1][2] = a * (-qx*qy); Q[1][3] = a * (-qx*qz);
        Q[2][0] = a * (-qw*qy);      Q[2][1] = a * (-qx*qy);      Q[2][2] = a * (1.0 - qy*qy); Q[2][3] = a * (-qy*qz);
        Q[3][0] = a * (-qw*qz);      Q[3][1] = a * (-qx*qz);      Q[3][2] = a * (-qy*qz);      Q[3][3] = a * (1.0 - qz*qz);

        // 9轴中 bz 也可观测 (磁力计约束yaw)，使用相同的零偏噪声
        double b = sigma_b_ * sigma_b_ * dt;
        Q[4][4] = b; Q[5][5] = b; Q[6][6] = b;

        // 协方差预测: P = F*P*F^T + Q
        double FP[7][7];
        mat_mul_7x7(F, P_, FP);
        mat_mul_7x7_AT(FP, F, P_);
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < 7; j++)
                P_[i][j] += Q[i][j];

        x_[0] = qwp; x_[1] = qxp; x_[2] = qyp; x_[3] = qzp;
        x_[4] = bxp; x_[5] = byp; x_[6] = bzp;
    }

    // ===================== 加速度计更新 =====================
    void ekf_update_accel(double ax, double ay, double az)
    {
        if (!ekf_inited_) return;

        double an2 = ax * ax + ay * ay + az * az;
        double an  = std::sqrt(an2);
        if (an < 0.5) return;
        if (accel_gate_ > 0.0 && std::fabs(an - expected_g_norm_) > accel_gate_) return;

        double qw = x_[0], qx = x_[1], qy = x_[2], qz = x_[3];
        double axn = ax / an, ayn = ay / an, azn = az / an;

        // 预测量测 h(q) — 重力方向在机体系的投影
        double hx = 2.0 * (qx * qz - qw * qy);
        double hy = 2.0 * (qw * qx + qy * qz);
        double hz = qw * qw - qx * qx - qy * qy + qz * qz;

        double y0 = axn - hx, y1 = ayn - hy, y2 = azn - hz;

        // 量测雅可比 H (3x7)
        double H[3][7];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 7; j++)
                H[i][j] = 0.0;
        H[0][0] = -2.0 * qy;  H[0][1] =  2.0 * qz;  H[0][2] = -2.0 * qw;  H[0][3] =  2.0 * qx;
        H[1][0] =  2.0 * qx;  H[1][1] =  2.0 * qw;  H[1][2] =  2.0 * qz;  H[1][3] =  2.0 * qy;
        H[2][0] =  2.0 * qw;  H[2][1] = -2.0 * qx;  H[2][2] = -2.0 * qy;  H[2][3] =  2.0 * qz;

        // R = σ_a² * I₃
        double ra = sigma_a_ * sigma_a_;
        double R[3][3] = {{ra, 0, 0}, {0, ra, 0}, {0, 0, ra}};

        kalman_update_3d(H, R, y0, y1, y2);
    }

    // ===================== 磁力计更新 =====================
    void ekf_update_mag(double mx, double my, double mz)
    {
        if (!ekf_inited_) return;

        double mn = std::sqrt(mx * mx + my * my + mz * mz);
        if (mn < 0.5) return;

        // 磁干扰检测: 幅值偏离参考过大则跳过
        double mag_ref_norm = std::sqrt(mag_ref_[0]*mag_ref_[0] +
                                        mag_ref_[1]*mag_ref_[1] +
                                        mag_ref_[2]*mag_ref_[2]);
        if (mag_gate_ > 0.0 && std::fabs(mn - mag_ref_norm) > mag_gate_ * mag_ref_norm) {
            return;
        }

        double qw = x_[0], qx = x_[1], qy = x_[2], qz = x_[3];
        double mxn = mx / mn, myn = my / mn, mzn = mz / mn;

        double nx = mag_ref_[0], ny = mag_ref_[1], nz = mag_ref_[2];

        // 预测量测 h_mag(q) = R(q)^T * mag_ref
        // R^T * n = body-frame 地磁预测

        double hx = (qw * qw + qx * qx - qy * qy - qz * qz)*nx + 2.0*(qx*qy + qw*qz)*ny + 2.0*(qx*qz - qw*qy)*nz;
        double hy = 2.0*(qx*qy - qw*qz)*nx + (qw * qw - qx * qx + qy * qy - qz * qz)*ny + 2.0*(qy*qz + qw*qx)*nz;
        double hz = 2.0*(qx*qz + qw*qy)*nx + 2.0*(qy*qz - qw*qx)*ny + (qw * qw - qx * qx - qy * qy + qz * qz)*nz;
        double y0 = mxn - hx, y1 = myn - hy, y2 = mzn - hz;

        // 量测雅可比 H_mag (3x7): ∂h_mag/∂x
        double H[3][7];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 7; j++)
                H[i][j] = 0.0;

        // ∂hx/∂q
        H[0][0] =  2.0*( qz*ny - qy*nz);
        H[0][1] =  2.0*( qy*ny + qz*nz);
        H[0][2] =  2.0*(-2.0*qy*nx + qx*ny - qw*nz);
        H[0][3] =  2.0*(-2.0*qz*nx + qw*ny + qx*nz);

        // ∂hy/∂q
        H[1][0] =  2.0*(-qz*nx + qx*nz);
        H[1][1] =  2.0*( qy*nx - 2.0*qx*ny + qw*nz);
        H[1][2] =  2.0*( qx*nx + qz*nz);
        H[1][3] =  2.0*(-qw*nx - 2.0*qz*ny + qy*nz);

        // ∂hz/∂q
        H[2][0] =  2.0*( qy*nx - qx*ny);
        H[2][1] =  2.0*( qz*nx - qw*ny - 2.0*qx*nz);
        H[2][2] =  2.0*( qw*nx + qz*ny - 2.0*qy*nz);
        H[2][3] =  2.0*( qx*nx + qy*ny);

        // R_mag = σ_m² * I₃
        double rm = sigma_m_ * sigma_m_;
        double R[3][3] = {{rm, 0, 0}, {0, rm, 0}, {0, 0, rm}};

        kalman_update_3d(H, R, y0, y1, y2);
    }

    // ===================== 通用 3D Kalman 更新 =====================
    void kalman_update_3d(const double H[3][7], const double R[3][3],
                          double y0, double y1, double y2)
    {
        double qw = x_[0], qx = x_[1], qy = x_[2], qz = x_[3];

        // S = H*P*H^T + R
        double HP[3][7];
        mat_mul_3x7_7x7(H, P_, HP);

        double S[3][3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                double sum = 0.0;
                for (int k = 0; k < 7; k++) sum += HP[i][k] * H[j][k];
                S[i][j] = sum + R[i][j];
            }

        double invS[3][3];
        if (!mat_inv_3x3(S, invS)) return;

        // K = P*H^T * S^{-1} = HP^T * S^{-1} (7x3)
        double K[7][3];
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < 3; j++) {
                double sum = 0.0;
                for (int k = 0; k < 3; k++) sum += HP[k][i] * invS[k][j];
                K[i][j] = sum;
            }

        // 状态更新
        x_[0] = qw + K[0][0]*y0 + K[0][1]*y1 + K[0][2]*y2;
        x_[1] = qx + K[1][0]*y0 + K[1][1]*y1 + K[1][2]*y2;
        x_[2] = qy + K[2][0]*y0 + K[2][1]*y1 + K[2][2]*y2;
        x_[3] = qz + K[3][0]*y0 + K[3][1]*y1 + K[3][2]*y2;
        x_[4] += K[4][0]*y0 + K[4][1]*y1 + K[4][2]*y2;
        x_[5] += K[5][0]*y0 + K[5][1]*y1 + K[5][2]*y2;
        x_[6] += K[6][0]*y0 + K[6][1]*y1 + K[6][2]*y2;

        // 四元数归一化
        double qn = inv_sqrt(x_[0]*x_[0] + x_[1]*x_[1] + x_[2]*x_[2] + x_[3]*x_[3]);
        x_[0] *= qn; x_[1] *= qn; x_[2] *= qn; x_[3] *= qn;

        // 协方差更新 (Joseph form)
        double IKH[7][7];
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < 7; j++) {
                double s = 0.0;
                for (int k = 0; k < 3; k++) s += K[i][k] * H[k][j];
                IKH[i][j] = (i == j ? 1.0 : 0.0) - s;
            }

        double P_temp[7][7];
        mat_mul_7x7(IKH, P_, P_temp);
        mat_mul_7x7_AT(P_temp, IKH, P_);

        // P += K*R*K^T
        double KR[7][3];
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < 3; j++) {
                double sum = 0.0;
                for (int k = 0; k < 3; k++) sum += K[i][k] * R[k][j];
                KR[i][j] = sum;
            }
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < 7; j++) {
                double sum = 0.0;
                for (int k = 0; k < 3; k++) sum += KR[i][k] * K[j][k];
                P_[i][j] += sum;
            }
    }

    // ===================== 四元数 → 欧拉角 (ZXY) =====================
    static void quat_to_euler(double w, double x, double y, double z,
                              double& roll, double& pitch, double& yaw)
    {
        double sinp = 2.0 * (w * x + y * z);
        if (sinp > 1.0)  sinp = 1.0;
        if (sinp < -1.0) sinp = -1.0;
        pitch = std::asin(sinp);

        roll = -std::atan2(2.0 * (x * z - w * y),
                           w * w - x * x - y * y + z * z);

        yaw = -std::atan2(2.0 * (x * y - w * z),
                          w * w - x * x + y * y - z * z);
    }

    // ===================== 发布与广播 =====================
    void publish_attitude(rclcpp::Time stamp)
    {
        geometry_msgs::msg::Vector3Stamped euler_msg;
        euler_msg.header.stamp = stamp;
        euler_msg.header.frame_id = "imu_link";

        double roll, pitch, yaw;
        quat_to_euler(x_[0], x_[1], x_[2], x_[3], roll, pitch, yaw);
        euler_msg.vector.x = roll;
        euler_msg.vector.y = pitch;
        euler_msg.vector.z = yaw;
        pub_attitude_->publish(euler_msg);

        static int cnt = 0;
        if (++cnt % 100 == 0) {
            RCLCPP_INFO(this->get_logger(),
                "Roll:%6.2f° Pitch:%6.2f° Yaw:%6.2f° | bias:[%+.4f %+.4f %+.4f]",
                roll * RAD_TO_DEG, pitch * RAD_TO_DEG, yaw * RAD_TO_DEG,
                x_[4], x_[5], x_[6]);
        }
    }

    void broadcast_tf(rclcpp::Time stamp)
    {
        geometry_msgs::msg::TransformStamped tf;
        tf.header.stamp = stamp;
        tf.header.frame_id = "world";
        tf.child_frame_id = "imu_attitude";
        tf.transform.rotation.w = x_[0];
        tf.transform.rotation.x = x_[1];
        tf.transform.rotation.y = x_[2];
        tf.transform.rotation.z = x_[3];
        tf_broadcaster_->sendTransform(tf);
    }

    // ---- 参数 ----
    std::string imu_topic_, mag_topic_, euler_topic_, output_topic_;
    double init_static_time_;
    double sigma_g_, sigma_b_, sigma_a_, sigma_m_;
    double accel_lpf_freq_, gyro_lpf_freq_, mag_lpf_freq_;
    double expected_g_norm_, accel_gate_, mag_gate_;
    bool enable_lpf_;

    // ---- 低通滤波器 ----
    bool lpf_configured_;
    LowPassFilter lpf_ax_, lpf_ay_, lpf_az_;
    LowPassFilter lpf_gx_, lpf_gy_, lpf_gz_;
    LowPassFilter lpf_mx_, lpf_my_, lpf_mz_;

    // ---- EKF 状态 ----
    double x_[7];       // [qw, qx, qy, qz, bx, by, bz]
    double P_[7][7];    // 协方差
    bool ekf_inited_;

    // ---- 磁力计 ----
    double mag_ref_[3];      // 导航系地磁参考 [mx, my, mz]
    double raw_mx_, raw_my_, raw_mz_;
    bool mag_received_;

    // ---- 硬件欧拉角 (用于初始yaw) ----
    double hw_roll_ = 0.0, hw_pitch_ = 0.0, hw_yaw_ = 0.0;
    bool euler_received_ = false;

    // ---- 初始化 ----
    bool initialized_;
    int init_sample_count_;
    double init_ax_sum_, init_ay_sum_, init_az_sum_;
    double init_mx_sum_, init_my_sum_, init_mz_sum_;

    // ---- ROS ----
    rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr sub_imu_;
    rclcpp::Subscription<sensor_msgs::msg::MagneticField>::SharedPtr sub_mag_;
    rclcpp::Subscription<geometry_msgs::msg::Vector3Stamped>::SharedPtr sub_euler_;
    rclcpp::Publisher<geometry_msgs::msg::Vector3Stamped>::SharedPtr pub_attitude_;
    rclcpp::Publisher<sensor_msgs::msg::Imu>::SharedPtr pub_filtered_imu_;
    std::shared_ptr<tf2_ros::TransformBroadcaster> tf_broadcaster_;
    rclcpp::Time last_time_;
};

} // namespace imu_9dof

int main(int argc, char** argv)
{
    rclcpp::init(argc, argv);
    auto node = std::make_shared<imu_9dof::EKF9DOF>();
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}
