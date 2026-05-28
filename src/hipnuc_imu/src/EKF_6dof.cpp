/*======================================================================
 *  EKF 6轴 AHRS (无磁力计)
 *
 *  导航系:  东北天 (ENU)      X=East, Y=North, Z=Up
 *  机体系:  右前上 (RFU)      X=Right, Y=Forward, Z=Up
 *  旋转:    ZXY 内旋          yaw(绕Z) → pitch(绕X) → roll(绕Y)
 *  四元数:  q_0=w, q_1=x, q_2=y, q_3=z (Hamilton)
 *  欧拉角:  roll=绕Y, pitch=绕X, yaw=绕Z
 *
 *  状态向量 (7D): [qw, qx, qy, qz, bias_x, bias_y, bias_z]
 *  预测量: 陀螺仪(角速度)驱动四元数运动学
 *  观测量: 加速度计(重力方向)修正roll和pitch
 *
 *  注意: 6轴无磁力计时, yaw及z轴陀螺零偏不可观测, 会随时间漂移
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
 * 给定四元数 [q_0, q_1, q_2, q_3]，方向余弦矩阵为：
C_{b2n} =
[q_0^2 + q_1^2 - q_2^2 - q_3^2 & 2(q_1q_2 - q_0q_3) & 2(q_1q_3 + q_0q_2) ;
2(q_1q_2 + q_0q_3) & q_0^2 - q_1^2 + q_2^2 - q_3^2 & 2(q_2q_3 - q_0q_1) ;
2(q_1q_3 - q_0q_2) & 2(q_2q_3 + q_0q_1) & q_0^2 - q_1^2 - q_2^2 + q_3^2]

 *======================================================================*/
//================================================================================
//  EKF 6轴 AHRS 姿态解算节点 (无磁力计)
//
//  导航系:  东北天 (ENU)      X=East, Y=North, Z=Up
//  机体系:  右前上 (RFU)      X=Right, Y=Forward, Z=Up
//  旋转:    ZXY 内旋          yaw(绕Z) → pitch(绕X) → roll(绕Y)
//  四元数:  q_0=w, q_1=x, q_2=y, q_3=z (Hamilton)
//  欧拉角:  roll=绕Y, pitch=绕X, yaw=绕Z
//
//  状态向量 (7D): [qw, qx, qy, qz, bias_x, bias_y, bias_z]
//  预测量: 陀螺仪(角速度)驱动四元数运动学
//  观测量: 加速度计(重力方向)修正roll和pitch
//
//  ⚠️ 注意: 6轴无磁力计时, yaw及z轴陀螺零偏不可观测, 会随时间漂移
//
//  ==================== 架构说明 ====================
//  ● ROS 接口层：参数声明、话题订阅/发布、TF广播、时间管理
//  ● 滤波器算法层：低通滤波器(EMA)、EKF预测/更新、矩阵运算、四元数/欧拉角转换
//  ● 移植时只需替换 ROS 接口（话题IO、参数读取、时间获取、日志打印）即可复用核心算法
//================================================================================

#include <memory>
#include <cmath>
#include <chrono>
#include <cstdint>

#include "rclcpp/rclcpp.hpp"
#include "sensor_msgs/msg/imu.hpp"
#include "geometry_msgs/msg/vector3_stamped.hpp"
#include "tf2_ros/transform_broadcaster.h"
#include "tf2_geometry_msgs/tf2_geometry_msgs.hpp"

using namespace std::chrono_literals;

/*======================================================================
 *  欧拉角/四元数定义及转换公式 (详见文件头顶)
 *======================================================================*/

// ========== 常量 ==========
#define DEG_TO_RAD  0.017453292519943295769236907684886
#define RAD_TO_DEG  57.295779513082320876798154814105

namespace imu_6dof
{

// ==================== 滤波器算法层 ====================
// 以下函数和类为纯数学/信号处理算法，与ROS无关，可直接移植

// ---------- 快速开方倒数 (用于四元数归一化) ----------
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

// ---------- 一阶低通滤波器 (指数移动平均) ----------
// 移植时无需修改，仅依赖 cutoff_freq 和 dt
class LowPassFilter
{
public:
    LowPassFilter() : initialized_(false), value_(0.0) {}

    // 设置截止频率和采样间隔
    void set_params(double cutoff_freq, double dt)
    {
        if (cutoff_freq <= 0.0) {
            alpha_ = 1.0;                     // 无滤波
        } else {
            double rc = 1.0 / (2.0 * M_PI * cutoff_freq);
            alpha_ = dt / (dt + rc);          // 平滑系数 = dt / (dt + RC)
            alpha_ = std::max(0.01, std::min(1.0, alpha_));
        }
    }

    // 更新滤波器，返回当前滤波值
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

// ---------- 小矩阵运算辅助 (7x7, 3x7, 3x3 求逆) ----------
// 用于EKF协方差计算

// 7x7 矩阵乘法: C = A * B
static void mat_mul_7x7(const double A[7][7], const double B[7][7], double C[7][7])
{
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 7; j++) {
            double sum = 0.0;
            for (int k = 0; k < 7; k++) sum += A[i][k] * B[k][j];
            C[i][j] = sum;
        }
}

// 7x7 * 7x7^T: C = A * B^T
static void mat_mul_7x7_AT(const double A[7][7], const double B[7][7], double C[7][7])
{
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 7; j++) {
            double sum = 0.0;
            for (int k = 0; k < 7; k++) sum += A[i][k] * B[j][k];
            C[i][j] = sum;
        }
}

// 3x7 * 7x7: C = A * B
static void mat_mul_3x7_7x7(const double A[3][7], const double B[7][7], double C[3][7])
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 7; j++) {
            double sum = 0.0;
            for (int k = 0; k < 7; k++) sum += A[i][k] * B[k][j];
            C[i][j] = sum;
        }
}

// 3x3 矩阵求逆 (用于新息协方差求逆)
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


// ==================== ROS 接口层 ====================
// 以下 EKF6DOF 类封装了ROS节点，算法核心部分可单独抽取移植

class EKF6DOF : public rclcpp::Node
{
public:
    EKF6DOF() : Node("ekf_6dof_node")
    {
        // ---------- ROS 参数声明 ----------
        this->declare_parameter<std::string>("imu_topic", "/IMU_data");
        this->declare_parameter<std::string>("output_topic", "/attitude");
        this->declare_parameter<double>("init_static_time", 2.0);

        // EKF噪声参数 (算法核心)
        this->declare_parameter<double>("sigma_g", 0.02);       // 陀螺仪角速度噪声 rad/s/√Hz
        this->declare_parameter<double>("sigma_b", 1.0e-4);     // 零偏随机游走 rad/s/√s/√Hz
        this->declare_parameter<double>("sigma_a", 0.01);       // 加速度计量测噪声 (归一化后)

        // 低通滤波参数 (算法核心)
        this->declare_parameter<double>("accel_lpf_freq", 30.0);
        this->declare_parameter<double>("gyro_lpf_freq", 30.0);
        this->declare_parameter<bool>("enable_lpf", false);      // 默认关闭，EKF自身有滤波效果

        // 加速度计量测预处理 (算法核心)
        this->declare_parameter<double>("expected_g_norm", 9.81);
        this->declare_parameter<double>("accel_gate", 2.0);

        // ---------- 读取参数 ----------
        this->get_parameter("imu_topic", imu_topic_);
        this->get_parameter("output_topic", output_topic_);
        this->get_parameter("init_static_time", init_static_time_);
        this->get_parameter("sigma_g", sigma_g_);
        this->get_parameter("sigma_b", sigma_b_);
        this->get_parameter("sigma_a", sigma_a_);
        this->get_parameter("accel_lpf_freq", accel_lpf_freq_);
        this->get_parameter("gyro_lpf_freq", gyro_lpf_freq_);
        this->get_parameter("enable_lpf", enable_lpf_);
        this->get_parameter("expected_g_norm", expected_g_norm_);
        this->get_parameter("accel_gate", accel_gate_);

        // ---------- ROS 通信接口 ----------
        sub_imu_ = this->create_subscription<sensor_msgs::msg::Imu>(
            imu_topic_, rclcpp::SensorDataQoS(),
            std::bind(&EKF6DOF::imu_callback, this, std::placeholders::_1));

        pub_attitude_ = this->create_publisher<geometry_msgs::msg::Vector3Stamped>(output_topic_, 10);
        pub_filtered_imu_ = this->create_publisher<sensor_msgs::msg::Imu>("/imu_filtered", 10);
        tf_broadcaster_ = std::make_shared<tf2_ros::TransformBroadcaster>(this);

        // ---------- 算法状态初始化 ----------
        // EKF 状态: 初始四元数为单位四元数，零偏为0
        x_[0] = 1.0; x_[1] = 0.0; x_[2] = 0.0; x_[3] = 0.0;
        x_[4] = 0.0; x_[5] = 0.0; x_[6] = 0.0;

        ekf_inited_ = false;
        initialized_ = false;
        init_sample_count_ = 0;
        init_ax_sum_ = init_ay_sum_ = init_az_sum_ = 0.0;
        lpf_configured_ = false;

        last_time_ = this->now();

        RCLCPP_INFO(this->get_logger(),
            "6-DOF EKF started | sigma_g=%.3f sigma_b=%.1e sigma_a=%.3f | "
            "g_norm=%.2f gate=%.2f | LPF %s (accel=%.0fHz, gyro=%.0fHz)",
            sigma_g_, sigma_b_, sigma_a_,
            expected_g_norm_, accel_gate_,
            enable_lpf_ ? "ON" : "OFF", accel_lpf_freq_, gyro_lpf_freq_);
    }

private:
    // ===================== ROS 回调：IMU 数据入口 =====================
    void imu_callback(const sensor_msgs::msg::Imu::SharedPtr msg)
    {
        // 获取 ROS 时间戳并计算 dt
        auto now = this->now();
        double dt = (now - last_time_).seconds();
        if (dt <= 0.0 || dt > 0.1) {
            last_time_ = now;
            return;
        }

        // ---------- 低通滤波器参数配置 (首次收到数据时，依赖于 dt) ----------
        if (!lpf_configured_ && enable_lpf_) {
            lpf_ax_.set_params(accel_lpf_freq_, dt);
            lpf_ay_.set_params(accel_lpf_freq_, dt);
            lpf_az_.set_params(accel_lpf_freq_, dt);
            lpf_gx_.set_params(gyro_lpf_freq_, dt);
            lpf_gy_.set_params(gyro_lpf_freq_, dt);
            lpf_gz_.set_params(gyro_lpf_freq_, dt);
            lpf_configured_ = true;
        }

        // 读取原始 IMU 数据 (RFU 坐标系)
        double gx = msg->angular_velocity.x;
        double gy = msg->angular_velocity.y;
        double gz = msg->angular_velocity.z;
        double ax = msg->linear_acceleration.x;
        double ay = msg->linear_acceleration.y;
        double az = msg->linear_acceleration.z;

        // ---------- 算法层：低通滤波 (可选) ----------
        if (enable_lpf_) {
            ax = lpf_ax_.update(ax);
            ay = lpf_ay_.update(ay);
            az = lpf_az_.update(az);
            gx = lpf_gx_.update(gx);
            gy = lpf_gy_.update(gy);
            gz = lpf_gz_.update(gz);
        }

        // ---------- ROS 发布：滤波后的 IMU 数据 ----------
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

        // ---------- 算法层：静态初始化 (利用静止时刻的加速度计算初始 roll/pitch) ----------
        if (!initialized_) {
            init_ax_sum_ += ax;
            init_ay_sum_ += ay;
            init_az_sum_ += az;
            init_sample_count_++;
            if (init_sample_count_ * dt >= init_static_time_) {
                double avg_ax = init_ax_sum_ / init_sample_count_;
                double avg_ay = init_ay_sum_ / init_sample_count_;
                double avg_az = init_az_sum_ / init_sample_count_;
                init_attitude(avg_ax, avg_ay, avg_az);   // 算法核心：由重力方向初始化姿态
                initialized_ = true;
                RCLCPP_INFO(this->get_logger(),
                    "EKF initialized (%d samples)", init_sample_count_);
            }
            last_time_ = now;
            return;
        }

        // ---------- 算法层：EKF 预测 & 更新 ----------
        ekf_predict(gx, gy, gz, dt);    // 状态预测 (四元数运动学 + 零偏随机游走)
        ekf_update(ax, ay, az);         // 状态更新 (加速度计观测重力方向)

        // ---------- ROS 输出：发布姿态角 + TF ----------
        publish_attitude(now);
        broadcast_tf(now);
        last_time_ = now;
    }

    // ===================== 算法层：初始对准 =====================
    // 输入: 静止时平均加速度 (m/s² 或 g 单位，归一化时使用)
    // 输出: 设置初始四元数和协方差矩阵
    void init_attitude(double ax, double ay, double az)
    {
        double norm = std::sqrt(ax * ax + ay * ay + az * az);
        if (norm < 1e-6) return;
        ax /= norm; ay /= norm; az /= norm;

        // 由重力方向计算 roll/pitch，yaw 设为 0
        double roll  = std::atan2(-ax, az);
        double pitch = std::atan2(ay, std::sqrt(ax * ax + az * az));
        double yaw = 0.0;

        RCLCPP_INFO(this->get_logger(),
            "Initial attitude: roll=%.2f°, pitch=%.2f°, yaw=%.2f°",
            roll * RAD_TO_DEG, pitch * RAD_TO_DEG, yaw * RAD_TO_DEG);

        // ZXY 欧拉角转四元数: q = qz(yaw) ⊗ qx(pitch) ⊗ qy(roll)
        double cy = std::cos(yaw * 0.5),   sy = std::sin(yaw * 0.5);
        double cp = std::cos(pitch * 0.5), sp = std::sin(pitch * 0.5);
        double cr = std::cos(roll * 0.5),  sr = std::sin(roll * 0.5);

        x_[0] = cy * cp * cr - sy * sp * sr;
        x_[1] = cy * sp * cr - sy * cp * sr;
        x_[2] = cy * cp * sr + sy * sp * cr;
        x_[3] = sy * cp * cr + cy * sp * sr;
        x_[4] = 0.0;  // bx
        x_[5] = 0.0;  // by
        x_[6] = 0.0;  // bz

        // 初始化协方差矩阵 (姿态部分 0.001，零偏部分 0.0001)
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < 7; j++)
                P_[i][j] = (i == j) ? ((i < 4) ? 0.001 : 0.0001) : 0.0;

        ekf_inited_ = true;
    }

    // ===================== 算法层：EKF 预测 =====================
    // 输入: 补偿前的角速度 (gx,gy,gz) 和采样间隔 dt
    // 更新: 状态 x_ 和协方差 P_
    void ekf_predict(double gx, double gy, double gz, double dt)
    {
        if (!ekf_inited_) return;

        double qw = x_[0], qx = x_[1], qy = x_[2], qz = x_[3];
        double bx = x_[4], by = x_[5], bz = x_[6];

        // 零偏补偿
        double wx = gx - bx;
        double wy = gy - by;
        double wz = gz - bz;

        // ---- 状态预测 (四元数运动学 Euler 积分) ----
        double half_dt = 0.5 * dt;
        double qwp = qw + half_dt * (-qx * wx - qy * wy - qz * wz);
        double qxp = qx + half_dt * ( qw * wx + qy * wz - qz * wy);
        double qyp = qy + half_dt * ( qw * wy - qx * wz + qz * wx);
        double qzp = qz + half_dt * ( qw * wz + qx * wy - qy * wx);

        // 归一化
        double qn = inv_sqrt(qwp * qwp + qxp * qxp + qyp * qyp + qzp * qzp);
        qwp *= qn; qxp *= qn; qyp *= qn; qzp *= qn;

        // 零偏保持不变
        double bxp = bx, byp = by, bzp = bz;

        // ---- 状态转移雅可比 F (7x7) ----
        double F[7][7];
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < 7; j++)
                F[i][j] = (i == j) ? 1.0 : 0.0;

        // ∂q̇/∂q (使用补偿后的角速度)
        F[0][1] = -half_dt * wx;  F[0][2] = -half_dt * wy;  F[0][3] = -half_dt * wz;
        F[1][0] =  half_dt * wx;  F[1][2] =  half_dt * wz;  F[1][3] = -half_dt * wy;
        F[2][0] =  half_dt * wy;  F[2][1] = -half_dt * wz;  F[2][3] =  half_dt * wx;
        F[3][0] =  half_dt * wz;  F[3][1] =  half_dt * wy;  F[3][2] = -half_dt * wx;

        // ∂q̇/∂b (使用未补偿前的四元数)
        F[0][4] =  half_dt * qx;  F[0][5] =  half_dt * qy;  F[0][6] =  half_dt * qz;
        F[1][4] = -half_dt * qw;  F[1][5] =  half_dt * qz;  F[1][6] = -half_dt * qy;
        F[2][4] = -half_dt * qz;  F[2][5] = -half_dt * qw;  F[2][6] =  half_dt * qx;
        F[3][4] =  half_dt * qy;  F[3][5] = -half_dt * qx;  F[3][6] = -half_dt * qw;

        // ---- 过程噪声协方差 Q (7x7) ----
        double Q[7][7];
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < 7; j++)
                Q[i][j] = 0.0;

        // 四元数部分: Q_q = (σ_g² * dt² / 4) * (I₄ - qqᵀ)
        double a = sigma_g_ * sigma_g_ * half_dt * half_dt;
        Q[0][0] = a * (1.0 - qw * qw);  Q[0][1] = a * (-qw * qx);
        Q[0][2] = a * (-qw * qy);       Q[0][3] = a * (-qw * qz);
        Q[1][0] = a * (-qw * qx);       Q[1][1] = a * (1.0 - qx * qx);
        Q[1][2] = a * (-qx * qy);       Q[1][3] = a * (-qx * qz);
        Q[2][0] = a * (-qw * qy);       Q[2][1] = a * (-qx * qy);
        Q[2][2] = a * (1.0 - qy * qy);  Q[2][3] = a * (-qy * qz);
        Q[3][0] = a * (-qw * qz);       Q[3][1] = a * (-qx * qz);
        Q[3][2] = a * (-qy * qz);       Q[3][3] = a * (1.0 - qz * qz);

        // 零偏随机游走: Q_b = σ_b² * dt (注意 bz 噪声极小，避免 yaw 漂移过快)
        double b_xy = sigma_b_ * sigma_b_ * dt;
        double b_z  = 1.0e-12 * dt;
        Q[4][4] = b_xy; Q[5][5] = b_xy; Q[6][6] = b_z;

        // ---- 协方差预测: P = F * P * Fᵀ + Q ----
        double FP[7][7];
        mat_mul_7x7(F, P_, FP);
        mat_mul_7x7_AT(FP, F, P_);
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < 7; j++)
                P_[i][j] += Q[i][j];

        // 保存预测状态
        x_[0] = qwp; x_[1] = qxp; x_[2] = qyp; x_[3] = qzp;
        x_[4] = bxp; x_[5] = byp; x_[6] = bzp;
    }

    // ===================== 算法层：EKF 更新 (加速度计观测) =====================
    // 输入: 加速度计测量值 (ax,ay,az)
    // 更新: 状态 x_ 和协方差 P_
    void ekf_update(double ax, double ay, double az)
    {
        if (!ekf_inited_) return;

        double qw = x_[0], qx = x_[1], qy = x_[2], qz = x_[3];

        // 加速度幅值检测 (拒收异常或大外部加速度)
        double an2 = ax * ax + ay * ay + az * az;
        double an  = std::sqrt(an2);
        if (an < 0.5) return;   // 幅值过小

        if (accel_gate_ > 0.0 && std::fabs(an - expected_g_norm_) > accel_gate_) {
            return;   // 偏离期望重力过大，存在明显线加速度
        }

        // 归一化加速度 (作为重力方向单位向量)
        double axn = ax / an, ayn = ay / an, azn = az / an;

        // ---- 观测模型 h(q) = 重力在机体系下的投影 ----
        double hx = 2.0 * (qx * qz - qw * qy);
        double hy = 2.0 * (qw * qx + qy * qz);
        double hz = qw * qw - qx * qx - qy * qy + qz * qz;

        // 新息
        double y0 = axn - hx;
        double y1 = ayn - hy;
        double y2 = azn - hz;

        // ---- 量测雅可比 H (3x7) ----
        double H[3][7];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 7; j++)
                H[i][j] = 0.0;

        H[0][0] = -2.0 * qy;  H[0][1] =  2.0 * qz;  H[0][2] = -2.0 * qw;  H[0][3] =  2.0 * qx;
        H[1][0] =  2.0 * qx;  H[1][1] =  2.0 * qw;  H[1][2] =  2.0 * qz;  H[1][3] =  2.0 * qy;
        H[2][0] =  2.0 * qw;  H[2][1] = -2.0 * qx;  H[2][2] = -2.0 * qy;  H[2][3] =  2.0 * qz;

        // ---- 量测噪声 R = σ_a² * I₃ ----
        double R[3][3] = {
            {sigma_a_ * sigma_a_, 0.0, 0.0},
            {0.0, sigma_a_ * sigma_a_, 0.0},
            {0.0, 0.0, sigma_a_ * sigma_a_}
        };

        // ---- S = H * P * Hᵀ + R ----
        double HP[3][7];
        mat_mul_3x7_7x7(H, P_, HP);

        double S[3][3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                double sum = 0.0;
                for (int k = 0; k < 7; k++) sum += HP[i][k] * H[j][k];
                S[i][j] = sum + R[i][j];
            }

        // 计算 S⁻¹
        double invS[3][3];
        if (!mat_inv_3x3(S, invS)) return;

        // ---- 卡尔曼增益 K = P * Hᵀ * S⁻¹ ----
        double K[7][3];
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < 3; j++) {
                double sum = 0.0;
                for (int k = 0; k < 3; k++) sum += HP[k][i] * invS[k][j];
                K[i][j] = sum;
            }

        // ---- 状态更新: x = x + K * y ----
        double dqw = K[0][0] * y0 + K[0][1] * y1 + K[0][2] * y2;
        double dqx = K[1][0] * y0 + K[1][1] * y1 + K[1][2] * y2;
        double dqy = K[2][0] * y0 + K[2][1] * y1 + K[2][2] * y2;
        double dqz = K[3][0] * y0 + K[3][1] * y1 + K[3][2] * y2;

        x_[0] = qw + dqw;
        x_[1] = qx + dqx;
        x_[2] = qy + dqy;
        x_[3] = qz + dqz;
        x_[4] += K[4][0] * y0 + K[4][1] * y1 + K[4][2] * y2;
        x_[5] += K[5][0] * y0 + K[5][1] * y1 + K[5][2] * y2;
        x_[6] += K[6][0] * y0 + K[6][1] * y1 + K[6][2] * y2;

        // 四元数归一化
        double qn = inv_sqrt(x_[0] * x_[0] + x_[1] * x_[1] +
                             x_[2] * x_[2] + x_[3] * x_[3]);
        x_[0] *= qn; x_[1] *= qn; x_[2] *= qn; x_[3] *= qn;

        // ---- 协方差更新 (Joseph form) ----
        // IKH = I - K*H
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

        // 加上 K*R*Kᵀ
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

    // ===================== 算法层：四元数转欧拉角 (ZXY 序列) =====================
    static void quat_to_euler(double w, double x, double y, double z,
                              double& roll, double& pitch, double& yaw)
    {
        double sinp = 2.0 * (w * x + y * z);
        if (sinp > 1.0)  sinp = 1.0;
        if (sinp < -1.0) sinp = -1.0;
        pitch = std::asin(sinp);

        roll  = -std::atan2(2.0 * (x * z - w * y),
                            w * w - x * x - y * y + z * z);
        yaw   = -std::atan2(2.0 * (x * y - w * z),
                            w * w - x * x + y * y - z * z);
    }

    // ===================== ROS 发布：欧拉角消息 =====================
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
                "Roll: %6.2f°, Pitch: %6.2f°, Yaw: %6.2f° | bias: [%+.4f %+.4f %+.4f]",
                roll * RAD_TO_DEG, pitch * RAD_TO_DEG, yaw * RAD_TO_DEG,
                x_[4], x_[5], x_[6]);
        }
    }

    // ===================== ROS 广播：tf 变换 =====================
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

    // ---------------- 以下为参数与状态变量 (混合) ----------------
    // ROS 参数
    std::string imu_topic_, output_topic_;
    double init_static_time_;
    double sigma_g_, sigma_b_, sigma_a_;
    double accel_lpf_freq_, gyro_lpf_freq_;
    double expected_g_norm_, accel_gate_;
    bool enable_lpf_;

    // 低通滤波器实例 (算法层)
    bool lpf_configured_;
    LowPassFilter lpf_ax_, lpf_ay_, lpf_az_;
    LowPassFilter lpf_gx_, lpf_gy_, lpf_gz_;

    // EKF 状态与协方差 (算法层)
    double x_[7];       // [qw, qx, qy, qz, bx, by, bz]
    double P_[7][7];
    bool ekf_inited_;

    // 初始化相关
    bool initialized_;
    int init_sample_count_;
    double init_ax_sum_, init_ay_sum_, init_az_sum_;

    // ROS 通信句柄
    rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr sub_imu_;
    rclcpp::Publisher<geometry_msgs::msg::Vector3Stamped>::SharedPtr pub_attitude_;
    rclcpp::Publisher<sensor_msgs::msg::Imu>::SharedPtr pub_filtered_imu_;
    std::shared_ptr<tf2_ros::TransformBroadcaster> tf_broadcaster_;
    rclcpp::Time last_time_;
};

} // namespace imu_6dof

// ===================== 主函数 (ROS 节点入口) =====================
int main(int argc, char** argv)
{
    rclcpp::init(argc, argv);
    auto node = std::make_shared<imu_6dof::EKF6DOF>();
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}