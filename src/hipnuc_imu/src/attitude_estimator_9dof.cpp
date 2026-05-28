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
 *  预测量: 陀螺仪驱动四元数运动学
 *  观测量: 加速度计(重力方向) + 磁力计(地磁方向)
 *
 *  欧拉角转四元数(ZXY): q = qz(yaw) ⊗ qx(pitch) ⊗ qy(roll)
    q0 = cp * cr * cy - sp * sr * sy;
    q1 = cr * cy * sp - cp * sr * sy;
    q2 = cp * cy * sr + cr * sp * sy;
    q3 = cp * cr * sy + sp * sr * cy;
 *
 *  四元数转欧拉角(ZXY):
    pitch = arcsin(2*(w*x + y*z))
    roll  = -atan2(2*(x*z - w*y), w^2 - x^2 - y^2 + z^2)
    yaw   = -atan2(2*(x*y - w*z), w^2 - x^2 + y^2 - z^2)
 *  机体到导航余弦矩阵:
 * C_b2n = [ q0^2 + q1^2 - q2^2 - q3^2,   2*(q1*q2 - q0*q3),       2*(q1*q3 + q0*q2);
    2*(q1*q2 + q0*q3),           q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 - q0*q1);
    2*(q1*q3 - q0*q2),           2*(q2*q3 + q0*q1),       q0^2 - q1^2 - q2^2 + q3^2 ];
    导航到机体旋转矩阵:
    C_n2B =
 [cos(roll)*cos(yaw) - sin(pitch)*sin(roll)*sin(yaw), cos(roll)*sin(yaw) + cos(yaw)*sin(pitch)*sin(roll), -cos(pitch)*sin(roll)]
[                              -cos(pitch)*sin(yaw),                                cos(pitch)*cos(yaw),            sin(pitch)]
[cos(yaw)*sin(roll) + cos(roll)*sin(pitch)*sin(yaw), sin(roll)*sin(yaw) - cos(roll)*cos(yaw)*sin(pitch),  cos(pitch)*cos(roll)]
                                                                                                                                                                                                                                                                                               -cos(pitch)*sin(roll),           sin(pitch),                                                                                                                                                                                                                                                                                                  cos(pitch)*cos(roll)]

    *======================================================================*/
#include <memory>
#include <cmath>
#include <chrono>
#include <algorithm>

#include "rclcpp/rclcpp.hpp"
#include "sensor_msgs/msg/imu.hpp"
#include "sensor_msgs/msg/magnetic_field.hpp"
#include "geometry_msgs/msg/vector3_stamped.hpp"
#include "tf2_ros/transform_broadcaster.h"
#include "tf2_geometry_msgs/tf2_geometry_msgs.hpp"

using namespace std::chrono_literals;

namespace attitude_estimator
{

    class AttitudeEstimator : public rclcpp::Node
    {
    public:
        AttitudeEstimator() : Node("attitude_estimator")
        {
            this->declare_parameter<std::string>("imu_topic", "/IMU_data");
            this->declare_parameter<std::string>("mag_topic", "/magnetic_data");
            this->declare_parameter<std::string>("output_topic", "/attitude");
            this->declare_parameter<double>("kp", 0.8);
            this->declare_parameter<double>("kp_yaw", 1);
            this->declare_parameter<double>("ki", 0.005);
            this->declare_parameter<double>("beta", 0.01);
            this->declare_parameter<int>("filter_type", 0);
            this->declare_parameter<bool>("use_mag", true);

            this->get_parameter("imu_topic", imu_topic_);
            this->get_parameter("mag_topic", mag_topic_);
            this->get_parameter("output_topic", output_topic_);
            this->get_parameter("kp", kp_);
            this->get_parameter("kp_yaw", kp_yaw_);
            this->get_parameter("ki", ki_);
            this->get_parameter("beta", beta_);
            this->get_parameter("filter_type", filter_type_);
            this->get_parameter("use_mag", use_mag_);

            sub_imu_ = this->create_subscription<sensor_msgs::msg::Imu>(
                imu_topic_, rclcpp::SensorDataQoS(),
                std::bind(&AttitudeEstimator::imu_callback, this, std::placeholders::_1));

            if (use_mag_)
            {
                sub_mag_ = this->create_subscription<sensor_msgs::msg::MagneticField>(
                    mag_topic_, rclcpp::SensorDataQoS(),
                    std::bind(&AttitudeEstimator::mag_callback, this, std::placeholders::_1));
            }

            pub_attitude_ = this->create_publisher<geometry_msgs::msg::Vector3Stamped>(output_topic_, 10);
            tf_broadcaster_ = std::make_shared<tf2_ros::TransformBroadcaster>(this);

            q0_ = 1.0;
            q1_ = 0.0;
            q2_ = 0.0;
            q3_ = 0.0;
            exInt_ = eyInt_ = 0.0;
            yaw_int_ = 0.0;
            mag_received_ = false;
            attitude_initialized_ = false;
            last_time_ = this->now();

            RCLCPP_INFO(this->get_logger(),
                        "Started (kp:%.2f kp_yaw:%.2f ki:%.4f mag:%d)",
                        kp_, kp_yaw_, ki_, use_mag_);
        }

    private:
        void mag_callback(const sensor_msgs::msg::MagneticField::SharedPtr msg)
        {
            mag_x_ = msg->magnetic_field.x;
            mag_y_ = msg->magnetic_field.y;
            mag_z_ = msg->magnetic_field.z;
            mag_received_ = true;
        }

        void imu_callback(const sensor_msgs::msg::Imu::SharedPtr msg)
        {
            auto now = this->now();
            double dt = (now - last_time_).seconds();
            if (dt <= 0.0 || dt > 0.1)
            {
                last_time_ = now;
                return;
            }

            double gx = msg->angular_velocity.x;
            double gy = msg->angular_velocity.y;
            double gz = msg->angular_velocity.z;
            double ax = msg->linear_acceleration.x;
            double ay = msg->linear_acceleration.y;
            double az = msg->linear_acceleration.z;

            if (!attitude_initialized_)
            {
                if (use_mag_ && mag_received_)
                {
                    init_attitude(ax, ay, az, mag_x_, mag_y_, mag_z_);
                    attitude_initialized_ = true;
                    double r, p, y;
                    quat_to_euler(q0_, q1_, q2_, q3_, r, p, y);
                    RCLCPP_INFO(this->get_logger(),
                                "Init -> R:%.1f P:%.1f Y:%.1f",
                                r * 180 / M_PI, p * 180 / M_PI, y * 180 / M_PI);
                }
                else if (!use_mag_)
                {
                    init_attitude(ax, ay, az, 0, 0, 0);
                    attitude_initialized_ = true;
                }
                if (!attitude_initialized_)
                {
                    integrate_gyro_only(gx, gy, gz, dt);
                }
                publish_attitude(now);
                broadcast_tf(now);
                last_time_ = now;
                return;
            }

            mahony_update(gx, gy, gz, ax, ay, az, dt);
            // madgwick_update(gx, gy, gz, ax, ay, az, dt);
            publish_attitude(now);
            broadcast_tf(now);
            last_time_ = now;
        }

        // ===================== 初始化 =====================

        void init_attitude(double ax, double ay, double az,
                           double mx, double my, double mz)
        {
            double norm = sqrt(ax * ax + ay * ay + az * az);
            if (norm < 1e-6)
                {return; } 
            ax /= norm;
            ay /= norm;
            az /= norm;

            double roll = -atan2(ax, az);
            double pitch = atan2(ay, sqrt(ax * ax + az * az));
            double nm = sqrt(mx * mx + my * my + mz * mz);
            mx /= nm;
            my /= nm;
            mz /= nm;
            double cp = cos(pitch), sp = sin(pitch);
            double cr = cos(roll), sr = sin(roll);
            double Xh = cr * mx + sr * mz;
            double Yh = sp * sr * mx + cp * my - sp * cr * mz;
            double yaw = atan2(Xh, Yh);
            euler_to_quat(roll, pitch, yaw);
        }

        void euler_to_quat(double roll, double pitch, double yaw)
        {
            double cr = cos(roll * 0.5), sr = sin(roll * 0.5);
            double cp = cos(pitch * 0.5), sp = sin(pitch * 0.5);
            double cy = cos(yaw * 0.5), sy = sin(yaw * 0.5);
            q0_ = cy * cp * cr - sy * sp * sr;
            q1_ = cy * sp * cr - sy * cp * sr;
            q2_ = cy * cp * sr + sy * sp * cr;
            q3_ = sy * cp * cr + cy * sp * sr;
            normalize_quat();
        }

        void integrate_gyro_only(double gx, double gy, double gz, double dt)
        {
            double q0d = 0.5 * (-q1_ * gx - q2_ * gy - q3_ * gz);
            double q1d = 0.5 * (q0_ * gx + q2_ * gz - q3_ * gy);
            double q2d = 0.5 * (q0_ * gy - q1_ * gz + q3_ * gx);
            double q3d = 0.5 * (q0_ * gz + q1_ * gy - q2_ * gx);
            q0_ += q0d * dt;
            q1_ += q1d * dt;
            q2_ += q2d * dt;
            q3_ += q3d * dt;
            normalize_quat();
        }

        // ===================== 核心滤波器 =====================
        void mahony_update(double gx, double gy, double gz,
                           double ax, double ay, double az, double dt)
        {
            // ---- 1. 归一化加速度 ----
            double norm = sqrt(ax * ax + ay * ay + az * az);
            if (norm < 1e-6)
                return;
            ax /= norm;
            ay /= norm;
            az /= norm;

            // ---- 2. 加速度计修正 roll/pitch ----
            double vx = 2.0 * (q1_ * q3_ - q0_ * q2_);
            double vy = 2.0 * (q0_ * q1_ + q2_ * q3_);
            double vz = q0_ * q0_ - q1_ * q1_ - q2_ * q2_ + q3_ * q3_;

            double ex = ay * vz - az * vy;
            double ey = az * vx - ax * vz;

            exInt_ += ex * ki_ * dt;
            eyInt_ += ey * ki_ * dt;
            const double IL = 0.1;
            exInt_ = std::max(-IL, std::min(IL, exInt_));
            eyInt_ = std::max(-IL, std::min(IL, eyInt_));

            gx += kp_ * ex + exInt_;
            gy += kp_ * ey + eyInt_;

            // ---- 3. 磁力计修正 yaw → 修正 gz ----
            double mx = mag_x_, my = mag_y_, mz = mag_z_;
            double nm = sqrt(mx * mx + my * my + mz * mz);
            bool mag_ok = (use_mag_ && mag_received_ && nm > 1e-6);

            if (mag_ok)
            {
                mx /= nm;
                my /= nm;
                mz /= nm;

                // 从四元数提取 roll/pitch（用于倾斜补偿）
                double sinp = 2.0 * (q0_ * q1_ + q2_ * q3_);
                sinp = std::clamp(sinp, -1.0, 1.0);
                double pitch_est = asin(sinp);
                double roll_est = -atan2(2.0 * (q1_ * q3_ - q0_ * q2_),
                                         q0_ * q0_ - q1_ * q1_ - q2_ * q2_ + q3_ * q3_);

                double cp = cos(pitch_est), sp = sin(pitch_est);
                double cr = cos(roll_est), sr = sin(roll_est);

                // ZXY 倾斜补偿
                double Xh = cr * mx + sr * mz;
                double Yh = sp * sr * mx + cp * my - sp * cr * mz;

                double yaw_mag = atan2(Xh, Yh);
                double yaw_quat = -atan2(2.0 * (q1_ * q2_ - q0_ * q3_),
                                         q0_ * q0_ - q1_ * q1_ + q2_ * q2_ - q3_ * q3_);

                double yaw_err = yaw_mag - yaw_quat;
                yaw_err = atan2(sin(yaw_err), cos(yaw_err)); // wrap to [-π, π]

                // ★★★ 关键修正：与 roll/pitch 保持一致的结构 ★★★
                // 积分项估计 yaw 轴陀螺仪零偏（量纲：rad/s）
                yaw_int_ += yaw_err * ki_ * dt;
                yaw_int_ = std::max(-IL, std::min(IL, yaw_int_));

                // 比例项 + 积分项 → 直接修正 gz（量纲：rad/s）
                gz += kp_ * yaw_err + yaw_int_;
            }

            // ---- 4. 四元数积分（所有轴的修正已在 gyros 中体现） ----
            double q0d = 0.5 * (-q1_ * gx - q2_ * gy - q3_ * gz);
            double q1d = 0.5 * (q0_ * gx + q2_ * gz - q3_ * gy);
            double q2d = 0.5 * (q0_ * gy - q1_ * gz + q3_ * gx);
            double q3d = 0.5 * (q0_ * gz + q1_ * gy - q2_ * gx);
            q0_ += q0d * dt;
            q1_ += q1d * dt;
            q2_ += q2d * dt;
            q3_ += q3d * dt;
            normalize_quat();
        }

        // ===================== 工具函数 =====================

        inline void normalize_quat()
        {
            double n = sqrt(q0_ * q0_ + q1_ * q1_ + q2_ * q2_ + q3_ * q3_);
            if (n < 1e-6)
                return;
            q0_ /= n;
            q1_ /= n;
            q2_ /= n;
            q3_ /= n;
        }

        void quat_to_euler(double w, double x, double y, double z,
                           double &roll, double &pitch, double &yaw)
        {
            pitch = asin(std::clamp(2.0 * (w * x + y * z), -1.0, 1.0));
            roll = -atan2(2.0 * (x * z - w * y), w * w - x * x - y * y + z * z);
            yaw = -atan2(2.0 * (x * y - w * z), w * w - x * x + y * y - z * z);
        }

        void publish_attitude(rclcpp::Time stamp)
        {
            geometry_msgs::msg::Vector3Stamped msg;
            msg.header.stamp = stamp;
            msg.header.frame_id = "imu_link";
            double r, p, y;
            quat_to_euler(q0_, q1_, q2_, q3_, r, p, y);
            msg.vector.x = r;
            msg.vector.y = p;
            msg.vector.z = y;
            pub_attitude_->publish(msg);

            static int cnt = 0;
            if (++cnt % 50 == 0)
                RCLCPP_INFO(this->get_logger(),
                            "R:%6.2f° P:%6.2f° Y:%6.2f°",
                            r * 180 / M_PI, p * 180 / M_PI, y * 180 / M_PI);
        }

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

        std::string imu_topic_, mag_topic_, output_topic_;
        double kp_, kp_yaw_, ki_, beta_;
        int filter_type_;
        bool use_mag_;
        double mag_x_ = 0, mag_y_ = 0, mag_z_ = 0;
        bool mag_received_ = false;
        bool attitude_initialized_ = false;
        double q0_, q1_, q2_, q3_;
        double exInt_, eyInt_;
        double yaw_int_;
        rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr sub_imu_;
        rclcpp::Subscription<sensor_msgs::msg::MagneticField>::SharedPtr sub_mag_;
        rclcpp::Publisher<geometry_msgs::msg::Vector3Stamped>::SharedPtr pub_attitude_;
        std::shared_ptr<tf2_ros::TransformBroadcaster> tf_broadcaster_;
        rclcpp::Time last_time_;
    };

} // namespace attitude_estimator

int main(int argc, char **argv)
{
    rclcpp::init(argc, argv);
    rclcpp::spin(std::make_shared<attitude_estimator::AttitudeEstimator>());
    rclcpp::shutdown();
    return 0;
}