// src/imu_9dof_attitude.cpp
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
        this->declare_parameter<double>("zeta", 0.95);      // 漂移抑制
        this->declare_parameter<double>("sample_rate", 1000.0); // 采样率
        
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

        // EKF 参数
        this->declare_parameter<double>("ekf_q_noise", 0.001);  // 过程噪声标准差
        this->declare_parameter<double>("ekf_r_noise", 0.01);   // 观测噪声标准差
        double q_noise, r_noise;
        this->get_parameter("ekf_q_noise", q_noise);
        this->get_parameter("ekf_r_noise", r_noise);

        // 初始化状态（单位四元数）
        x_[0] = 1.0; x_[1] = 0.0; x_[2] = 0.0; x_[3] = 0.0;

        // 初始化协方差矩阵（单位矩阵）
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                P_[i][j] = (i == j) ? 0.1 : 0.0;  // 初始不确定性

        // 过程噪声矩阵（对角阵）
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                Q_[i][j] = (i == j) ? q_noise : 0.0;

        // 观测噪声矩阵（对角阵）
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                R_[i][j] = (i == j) ? r_noise : 0.0;

        ekf_initialized_ = true;
        RCLCPP_INFO(this->get_logger(), "EKF initialized with Q=%.4f, R=%.4f", q_noise, r_noise);
    }

private:
    bool bias_calibrated_;
    int bias_samples_;
    double gx_bias_, gy_bias_, gz_bias_;

    double x_[4];        // 四元数状态 [q0, q1, q2, q3]
    double P_[4][4];     // 状态协方差矩阵 4x4
    double Q_[4][4];     // 过程噪声协方差矩阵 4x4
    double R_[3][3];     // 观测噪声协方差矩阵 3x3
    bool ekf_initialized_;

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
        //complementary_filter(gx, gy, gz, ax, ay, az, dt);
        
        // ============================================
        // 方法2: Mahony 滤波 (更稳定，推荐)
        // 取消注释下面的代码，并注释掉 complementary_filter
        // ============================================
        //mahony_filter(gx, gy, gz, ax, ay, az, dt);
        
        // ============================================
        // 方法3: Madgwick 滤波 (计算量稍大)
        // ============================================
        // madgwick_filter(gx, gy, gz, ax, ay, az, dt);
        
        // ============================================
        // 方法4: EKF 滤波 (计算量稍大)
        // ============================================
        ekf_filter(gx, gy, gz, ax, ay, az, dt);
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
    
    void ekf_filter(double gx, double gy, double gz,
                double ax, double ay, double az,
                double dt)
    {
    if (!ekf_initialized_) return;

    // ========== 1. 预测步骤 ==========
    // 构建四元数更新矩阵 A = I + 0.5*dt*Ω(ω)
    double omega[4][4] = {0};
    omega[0][1] = -gx; omega[0][2] = -gy; omega[0][3] = -gz;
    omega[1][0] =  gx; omega[1][2] =  gz; omega[1][3] = -gy;
    omega[2][0] =  gy; omega[2][1] = -gz; omega[2][3] =  gx;
    omega[3][0] =  gz; omega[3][1] =  gy; omega[3][2] = -gx;

    double A[4][4];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            A[i][j] = (i == j ? 1.0 : 0.0) + 0.5 * dt * omega[i][j];
        }
    }

    // 状态预测：x_ = A * x_
    double new_x[4] = {0};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            new_x[i] += A[i][j] * x_[j];
    for (int i = 0; i < 4; ++i) x_[i] = new_x[i];

    // 协方差预测：P = A * P * A^T + Q
    double AP[4][4] = {0};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k)
                AP[i][j] += A[i][k] * P_[k][j];

    double APA[4][4] = {0};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k)
                APA[i][j] += AP[i][k] * A[j][k];  // A^T 的转置索引

    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            P_[i][j] = APA[i][j] + Q_[i][j];

    // 归一化四元数（防止数值漂移）
    double norm = sqrt(x_[0]*x_[0] + x_[1]*x_[1] + x_[2]*x_[2] + x_[3]*x_[3]);
    if (norm > 1e-6) {
        for (int i = 0; i < 4; ++i) x_[i] /= norm;
    }

    // ========== 2. 更新步骤（使用加速度计） ==========
    // 归一化加速度计测量值
    double norm_a = sqrt(ax*ax + ay*ay + az*az);
    if (norm_a < 1e-6) return;
    double z[3] = {ax / norm_a, ay / norm_a, az / norm_a};

    // 预测的重力向量（h(x)）
    double q0 = x_[0], q1 = x_[1], q2 = x_[2], q3 = x_[3];
    double h[3] = { 2*(q1*q3 - q0*q2),
                    2*(q0*q1 + q2*q3),
                    q0*q0 - q1*q1 - q2*q2 + q3*q3 };

    // 观测残差
    double y[3] = { z[0] - h[0], z[1] - h[1], z[2] - h[2] };

    // 观测矩阵 H (3x4) = ∂h/∂x
    double H[3][4] = {0};
    H[0][0] = -2*q2;   H[0][1] =  2*q3;   H[0][2] = -2*q0;   H[0][3] =  2*q1;
    H[1][0] =  2*q1;   H[1][1] =  2*q0;   H[1][2] =  2*q3;   H[1][3] =  2*q2;
    H[2][0] =  2*q0;   H[2][1] = -2*q1;   H[2][2] = -2*q2;   H[2][3] =  2*q3;

    // 卡尔曼增益 K = P * H^T * (H * P * H^T + R)^{-1}
    double PHt[4][3] = {0};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 4; ++k)
                PHt[i][j] += P_[i][k] * H[j][k];  // H 的转置是 H[j][k]

    double S[3][3] = {0};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double sum = 0;
            for (int k = 0; k < 4; ++k)
                sum += H[i][k] * PHt[k][j];
            S[i][j] = sum + R_[i][j];
        }

    // 求 S 的逆（3x3 对称正定矩阵）
    double invS[3][3];
    double det = S[0][0]* (S[1][1]*S[2][2] - S[2][1]*S[1][2])
               - S[0][1]* (S[1][0]*S[2][2] - S[1][2]*S[2][0])
               + S[0][2]* (S[1][0]*S[2][1] - S[1][1]*S[2][0]);
    if (fabs(det) < 1e-12) return;
    double invDet = 1.0 / det;
    invS[0][0] = (S[1][1]*S[2][2] - S[1][2]*S[2][1]) * invDet;
    invS[0][1] = (S[0][2]*S[2][1] - S[0][1]*S[2][2]) * invDet;
    invS[0][2] = (S[0][1]*S[1][2] - S[0][2]*S[1][1]) * invDet;
    invS[1][0] = (S[1][2]*S[2][0] - S[1][0]*S[2][2]) * invDet;
    invS[1][1] = (S[0][0]*S[2][2] - S[0][2]*S[2][0]) * invDet;
    invS[1][2] = (S[0][2]*S[1][0] - S[0][0]*S[1][2]) * invDet;
    invS[2][0] = (S[1][0]*S[2][1] - S[1][1]*S[2][0]) * invDet;
    invS[2][1] = (S[0][1]*S[2][0] - S[0][0]*S[2][1]) * invDet;
    invS[2][2] = (S[0][0]*S[1][1] - S[0][1]*S[1][0]) * invDet;

    double K[4][3];  // 卡尔曼增益 4x3
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j) {
            double sum = 0;
            for (int k = 0; k < 3; ++k)
                sum += PHt[i][k] * invS[k][j];
            K[i][j] = sum;
        }

    // 状态更新：x = x + K * y
    double delta[4] = {0};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            delta[i] += K[i][j] * y[j];
    for (int i = 0; i < 4; ++i) x_[i] += delta[i];

    // 归一化四元数
    norm = sqrt(x_[0]*x_[0] + x_[1]*x_[1] + x_[2]*x_[2] + x_[3]*x_[3]);
    if (norm > 1e-6) {
        for (int i = 0; i < 4; ++i) x_[i] /= norm;
    }

    // 协方差更新：P = (I - K*H) * P
    double KH[4][4] = {0};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 3; ++k)
                KH[i][j] += K[i][k] * H[k][j];

    double newP[4][4] = {0};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            newP[i][j] = (i == j ? 1.0 : 0.0) - KH[i][j];
        }
    double temp[4][4] = {0};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k)
                temp[i][j] += newP[i][k] * P_[k][j];
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            P_[i][j] = temp[i][j];

    // 将 EKF 估计的四元数同步到类全局四元数（用于发布及 TF）
    q0_ = x_[0];
    q1_ = x_[1];
    q2_ = x_[2];
    q3_ = x_[3];
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