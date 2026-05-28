// src/attitude_estimator.cpp
#include <memory>
#include <cmath>
#include <chrono>

#include "rclcpp/rclcpp.hpp"
#include "sensor_msgs/msg/imu.hpp"
#include "geometry_msgs/msg/quaternion.hpp"
#include "geometry_msgs/msg/vector3_stamped.hpp"
#include "tf2/LinearMath/Quaternion.h"
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
        // 参数声明
        this->declare_parameter<std::string>("imu_topic", "/IMU_data");
        this->declare_parameter<std::string>("output_topic", "/attitude");
        this->declare_parameter<double>("beta", 0.1);      // 滤波增益
        this->declare_parameter<double>("zeta", 0.8);      // 漂移抑制
        this->declare_parameter<double>("sample_rate", 200.0); // 采样率
        
        this->get_parameter("imu_topic", imu_topic_);
        this->get_parameter("output_topic", output_topic_);
        this->get_parameter("beta", beta_);
        this->get_parameter("zeta", zeta_);
        this->get_parameter("sample_rate", sample_rate_);
        
        // 订阅 IMU 数据
        subscription_ = this->create_subscription<sensor_msgs::msg::Imu>(
            imu_topic_, 
            rclcpp::SensorDataQoS(),
            std::bind(&AttitudeEstimator::imu_callback, this, std::placeholders::_1));
        
        // 发布姿态结果
        attitude_pub_ = this->create_publisher<geometry_msgs::msg::Vector3Stamped>(
            output_topic_, 10);
        
        // TF 广播器
        tf_broadcaster_ = std::make_shared<tf2_ros::TransformBroadcaster>(this);
        
        // 初始化四元数
        q0_ = 1.0;
        q1_ = 0.0;
        q2_ = 0.0;
        q3_ = 0.0;
        
        // 初始化陀螺仪漂移
        exInt_ = 0.0;
        eyInt_ = 0.0;
        ezInt_ = 0.0;
        
        last_time_ = this->now();
        
        RCLCPP_INFO(this->get_logger(), "Attitude Estimator initialized");
        RCLCPP_INFO(this->get_logger(), "Beta: %.3f, Zeta: %.3f", beta_, zeta_);
    }

private:
    void imu_callback(const sensor_msgs::msg::Imu::SharedPtr msg)
    {
        // 获取时间差（秒）
        auto current_time = this->now();
        double dt = (current_time - last_time_).seconds();
        if (dt <= 0 || dt > 0.1) {
            last_time_ = current_time;
            return;
        }
        
        // 提取数据
        double gx = msg->angular_velocity.x;   // rad/s
        double gy = msg->angular_velocity.y;
        double gz = msg->angular_velocity.z;
        
        double ax = msg->linear_acceleration.x; // m/s^2
        double ay = msg->linear_acceleration.y;
        double az = msg->linear_acceleration.z;
        
        // ============================================
        // 方法1: 互补滤波 (Complementary Filter)
        // ============================================
        // complementary_filter(gx, gy, gz, ax, ay, az, dt);
        
        // ============================================
        // 方法2: Mahony 滤波 (更稳定，推荐)
        // 取消注释下面的代码，并注释掉 complementary_filter
        // ============================================
        //mahony_filter(gx, gy, gz, ax, ay, az, dt);
        
        // ============================================
        // 方法3: Madgwick 滤波 (计算量稍大)
        // ============================================
        madgwick_filter(gx, gy, gz, ax, ay, az, dt);
        
        // 发布姿态
        publish_attitude(current_time);
        
        // 广播 TF
        broadcast_tf(current_time);
        
        last_time_ = current_time;
    }
    
    // 互补滤波算法
    void complementary_filter(double gx, double gy, double gz,
                              double ax, double ay, double az,
                              double dt)
    {
        // 归一化加速度计数据
        double norm = sqrt(ax*ax + ay*ay + az*az);
        if (norm < 0.001) return;
        
        ax /= norm;
        ay /= norm;
        az /= norm;
        
        // 从当前四元数计算重力向量
        double vx = 2*(q1_*q3_ - q0_*q2_);
        double vy = 2*(q0_*q1_ + q2_*q3_);
        double vz = q0_*q0_ - q1_*q1_ - q2_*q2_ + q3_*q3_;
        
        // 计算误差（加速度计和重力向量的叉积）
        double ex = ay*vz - az*vy;
        double ey = az*vx - ax*vz;
        double ez = ax*vy - ay*vx;
        
        // 积分误差补偿陀螺仪漂移
        exInt_ += ex * dt * beta_;
        eyInt_ += ey * dt * beta_;
        ezInt_ += ez * dt * beta_;
        
        // 修正陀螺仪数据
        gx += beta_ * ex + exInt_;
        gy += beta_ * ey + eyInt_;
        gz += beta_ * ez + ezInt_;
        
        // 四元数微分方程
        double q0_dot = (-q1_*gx - q2_*gy - q3_*gz) * 0.5;
        double q1_dot = ( q0_*gx - q3_*gy + q2_*gz) * 0.5;
        double q2_dot = ( q3_*gx + q0_*gy - q1_*gz) * 0.5;
        double q3_dot = (-q2_*gx + q1_*gy + q0_*gz) * 0.5;
        
        // 积分更新四元数
        q0_ += q0_dot * dt;
        q1_ += q1_dot * dt;
        q2_ += q2_dot * dt;
        q3_ += q3_dot * dt;
        
        // 归一化四元数
        norm = sqrt(q0_*q0_ + q1_*q1_ + q2_*q2_ + q3_*q3_);
        if (norm > 0.001) {
            q0_ /= norm;
            q1_ /= norm;
            q2_ /= norm;
            q3_ /= norm;
        }
    }
    
    // Mahony 滤波算法（更精确）
    void mahony_filter(double gx, double gy, double gz,
                       double ax, double ay, double az,
                       double dt)
    {
        double recipNorm;
        double vx, vy, vz;
        double ex, ey, ez;
        
        // 归一化加速度计
        recipNorm = 1.0 / sqrt(ax*ax + ay*ay + az*az);
        ax *= recipNorm;
        ay *= recipNorm;
        az *= recipNorm;
        
        // 估计重力方向
        vx = 2.0f * (q1_*q3_ - q0_*q2_);
        vy = 2.0f * (q0_*q1_ + q2_*q3_);
        vz = q0_*q0_ - q1_*q1_ - q2_*q2_ + q3_*q3_;
        
        // 误差 = 加速度计与重力向量的叉积
        ex = (ay*vz - az*vy);
        ey = (az*vx - ax*vz);
        ez = (ax*vy - ay*vx);
        
        // 积分误差
        exInt_ += ex * dt * zeta_;
        eyInt_ += ey * dt * zeta_;
        ezInt_ += ez * dt * zeta_;
        
        // 应用反馈
        gx += beta_ * ex + exInt_;
        gy += beta_ * ey + eyInt_;
        gz += beta_ * ez + ezInt_;
        
        // 四元数微分
        double q0_dot = (-q1_*gx - q2_*gy - q3_*gz) * 0.5f;
        double q1_dot = ( q0_*gx + q2_*gz - q3_*gy) * 0.5f;
        double q2_dot = ( q0_*gy - q1_*gz + q3_*gx) * 0.5f;
        double q3_dot = ( q0_*gz + q1_*gy - q2_*gx) * 0.5f;
        
        // 积分
        q0_ += q0_dot * dt;
        q1_ += q1_dot * dt;
        q2_ += q2_dot * dt;
        q3_ += q3_dot * dt;
        
        // 归一化
        recipNorm = 1.0 / sqrt(q0_*q0_ + q1_*q1_ + q2_*q2_ + q3_*q3_);
        q0_ *= recipNorm;
        q1_ *= recipNorm;
        q2_ *= recipNorm;
        q3_ *= recipNorm;
    }
    
    // Madgwick 滤波算法（高性能）
    void madgwick_filter(double gx, double gy, double gz,
                         double ax, double ay, double az,
                         double dt)
    {
        double q0 = q0_, q1 = q1_, q2 = q2_, q3 = q3_;
        double recipNorm;
        
        // 加速度计归一化
        recipNorm = 1.0 / sqrt(ax*ax + ay*ay + az*az);
        ax *= recipNorm;
        ay *= recipNorm;
        az *= recipNorm;
        
        // 梯度下降算法
        double f1, f2, f3, f4;
        double J_11or24, J_12or23, J_13or22, J_14or21;
        double J_32, J_33;
        
        f1 = 2*(q1*q3 - q0*q2) - ax;
        f2 = 2*(q0*q1 + q2*q3) - ay;
        f3 = 2*(0.5 - q1*q1 - q2*q2) - az;
        
        J_11or24 = 2*q2;
        J_12or23 = 2*q3;
        J_13or22 = 2*q0;
        J_14or21 = 2*q1;
        J_32 = 2*q1;
        J_33 = 2*q2;
        
        // 计算梯度
        double grad0 = J_14or21*f2 - J_11or24*f1;
        double grad1 = J_12or23*f1 + J_13or22*f2 + J_32*f3;
        double grad2 = J_12or23*f2 - J_14or21*f1 + J_33*f3;
        double grad3 = J_11or24*f2 - J_13or22*f1;
        
        // 归一化梯度
        recipNorm = 1.0 / sqrt(grad0*grad0 + grad1*grad1 + grad2*grad2 + grad3*grad3);
        grad0 *= recipNorm;
        grad1 *= recipNorm;
        grad2 *= recipNorm;
        grad3 *= recipNorm;
        
        // 陀螺仪积分和梯度下降的结合
        double qDot1 = 0.5 * (-q1*gx - q2*gy - q3*gz) - beta_ * grad0;
        double qDot2 = 0.5 * ( q0*gx + q2*gz - q3*gy) - beta_ * grad1;
        double qDot3 = 0.5 * ( q0*gy - q1*gz + q3*gx) - beta_ * grad2;
        double qDot4 = 0.5 * ( q0*gz + q1*gy - q2*gx) - beta_ * grad3;
        
        // 积分
        q0_ += qDot1 * dt;
        q1_ += qDot2 * dt;
        q2_ += qDot3 * dt;
        q3_ += qDot4 * dt;
        
        // 归一化
        recipNorm = 1.0 / sqrt(q0_*q0_ + q1_*q1_ + q2_*q2_ + q3_*q3_);
        q0_ *= recipNorm;
        q1_ *= recipNorm;
        q2_ *= recipNorm;
        q3_ *= recipNorm;
    }
    
    // 将四元数转换为欧拉角（Roll, Pitch, Yaw）
    void quaternion_to_euler(double q0, double q1, double q2, double q3,
                             double& roll, double& pitch, double& yaw)
    {
        // Roll (x-axis rotation)
        double sinr_cosp = 2.0 * (q0*q1 + q2*q3);
        double cosr_cosp = 1.0 - 2.0*(q1*q1 + q2*q2);
        roll = atan2(sinr_cosp, cosr_cosp);
        
        // Pitch (y-axis rotation)
        double sinp = 2.0 * (q0*q2 - q3*q1);
        if (fabs(sinp) >= 1.0)
            pitch = copysign(M_PI/2.0, sinp);
        else
            pitch = asin(sinp);
        
        // Yaw (z-axis rotation)
        double siny_cosp = 2.0 * (q0*q3 + q1*q2);
        double cosy_cosp = 1.0 - 2.0*(q2*q2 + q3*q3);
        yaw = atan2(siny_cosp, cosy_cosp);
    }
    
    void publish_attitude(rclcpp::Time stamp)
    {
        geometry_msgs::msg::Vector3Stamped euler_msg;
        euler_msg.header.stamp = stamp;
        euler_msg.header.frame_id = "imu_link";
        
        double roll, pitch, yaw;
        quaternion_to_euler(q0_, q1_, q2_, q3_, roll, pitch, yaw);
        
        euler_msg.vector.x = roll;   // 弧度
        euler_msg.vector.y = pitch;
        euler_msg.vector.z = yaw;
        
        attitude_pub_->publish(euler_msg);
        
        // 每100帧打印一次（避免刷屏）
        static int count = 0;
        if (++count % 100 == 0) {
            RCLCPP_INFO(this->get_logger(), 
                       "Roll: %.2f°, Pitch: %.2f°, Yaw: %.2f°",
                       roll * 180.0/M_PI, pitch * 180.0/M_PI, yaw * 180.0/M_PI);
        }
    }
    
    void broadcast_tf(rclcpp::Time stamp)
    {
        geometry_msgs::msg::TransformStamped transform;
        transform.header.stamp = stamp;
        transform.header.frame_id = "world";
        transform.child_frame_id = "imu_attitude";
        
        transform.transform.translation.x = 0.0;
        transform.transform.translation.y = 0.0;
        transform.transform.translation.z = 0.0;
        
        transform.transform.rotation.x = q1_;
        transform.transform.rotation.y = q2_;
        transform.transform.rotation.z = q3_;
        transform.transform.rotation.w = q0_;
        
        tf_broadcaster_->sendTransform(transform);
    }
    
    // 参数
    std::string imu_topic_;
    std::string output_topic_;
    double beta_;      // 滤波增益
    double zeta_;      // 漂移抑制
    double sample_rate_;
    
    // 四元数
    double q0_, q1_, q2_, q3_;
    
    // 误差积分
    double exInt_, eyInt_, ezInt_;
    
    // ROS 接口
    rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr subscription_;
    rclcpp::Publisher<geometry_msgs::msg::Vector3Stamped>::SharedPtr attitude_pub_;
    std::shared_ptr<tf2_ros::TransformBroadcaster> tf_broadcaster_;
    
    rclcpp::Time last_time_;
};

} // namespace attitude_estimator

int main(int argc, char** argv)
{
    rclcpp::init(argc, argv);
    auto node = std::make_shared<attitude_estimator::AttitudeEstimator>();
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}