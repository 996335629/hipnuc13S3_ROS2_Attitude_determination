#include <unistd.h>
#include <memory>
#include <iostream>
#include <iomanip>
#include "rclcpp/rclcpp.hpp"
#include "sensor_msgs/msg/imu.hpp"
#include "sensor_msgs/msg/magnetic_field.hpp"
#include "geometry_msgs/msg/vector3_stamped.hpp"

rclcpp::Node::SharedPtr nh = nullptr;
using namespace std;

// IMU 数据回调
void imu_callback(const sensor_msgs::msg::Imu::SharedPtr msg)
{
    static int count = 0;
    if (++count % 100 != 0)
        return;

    cout << "\n========== IMU Data ==========" << endl;
    cout << "header.stamp: sec=" << msg->header.stamp.sec 
         << " nsec=" << msg->header.stamp.nanosec << endl;
    //cout << "frame_id: " << msg->header.frame_id << endl;

    cout << "orientation (w,x,y,z): " 
         << fixed << setprecision(6)
         << msg->orientation.w << ", "
         << msg->orientation.x << ", "
         << msg->orientation.y << ", "
         << msg->orientation.z << endl;
    double w=msg->orientation.w;
    double x=msg->orientation.x;
    double y=msg->orientation.y;
    double z=msg->orientation.z;
    double pitch = asin(std::clamp(2.0 * (w * x + y * z), -1.0, 1.0))*57.29577951308;
    double roll = -atan2(2.0 * (x * z - w * y), w * w - x * x - y * y + z * z)*57.29577951308;
    double yaw = -atan2(2.0 * (x * y - w * z), w * w - x * x + y * y - z * z)*57.29577951308;

    cout << "quat_to_roll, pitch, yaw (°): " 
         << fixed << setprecision(6)
         << roll << ", "
         << pitch << ", "
         << yaw << endl;
    // cout << "angular_velocity (rad/s): "
    //      << msg->angular_velocity.x << ", "
    //      << msg->angular_velocity.y << ", "
    //      << msg->angular_velocity.z << endl;

    // cout << "linear_acceleration (m/s^2): "
    //      << msg->linear_acceleration.x << ", "
    //      << msg->linear_acceleration.y << ", "
    //      << msg->linear_acceleration.z << endl;
    cout << "================================" << endl;
}

// 欧拉角回调 (Vector3Stamped)
void euler_callback(const geometry_msgs::msg::Vector3Stamped::SharedPtr msg)
{
    // static int count = 0;
    // if (++count % 100 != 0)
    //     return;

    // cout << "\n========== Euler Data ==========" << endl;
    // cout << "stamp: sec=" << msg->header.stamp.sec 
    //      << " nsec=" << msg->header.stamp.nanosec << endl;
    // cout << "frame_id: " << msg->header.frame_id << endl;
    // cout << "roll, pitch, yaw (°): "
    //      << fixed << setprecision(6)
    //      << msg->vector.x*57.29577951308 << ", "
    //      << msg->vector.y*57.29577951308 << ", "
    //      << msg->vector.z*57.29577951308 << endl;
    // cout << "================================" << endl;
}

// 磁力计回调 (MagneticField)
void magnetic_callback(const sensor_msgs::msg::MagneticField::SharedPtr msg)
{
    static int count = 0;
    if (++count % 100 != 0)
        return;

    cout << "\n========== Magnetic Field ==========" << endl;
    cout << "stamp: sec=" << msg->header.stamp.sec 
         << " nsec=" << msg->header.stamp.nanosec << endl;
    //cout << "frame_id: " << msg->header.frame_id << endl;
    cout << "magnetic_field (tesla): "
         << fixed << setprecision(6)
         << msg->magnetic_field.x*1e6 << ", "
         << msg->magnetic_field.y*1e6 << ", "
         << msg->magnetic_field.z*1e6 << ", "
         << 1e6*sqrt(msg->magnetic_field.x*msg->magnetic_field.x+msg->magnetic_field.y*msg->magnetic_field.y+msg->magnetic_field.z*msg->magnetic_field.z) << endl;;
    // cout << "====================================" << endl;
}

// void quat_to_euler(double w, double x, double y, double z)
// {
//     pitch = asin(std::clamp(2.0 * (w * x + y * z), -1.0, 1.0));
//     roll = -atan2(2.0 * (x * z - w * y), w * w - x * x - y * y + z * z);
//     yaw = -atan2(2.0 * (x * y - w * z), w * w - x * x + y * y - z * z);
// }

int main(int argc, const char* argv[])
{
    rclcpp::init(argc, argv);
    nh = std::make_shared<rclcpp::Node>("imu_sub");

    // 使用 SensorDataQoS 保证数据实时性
    rclcpp::SensorDataQoS qos;

    // 订阅 IMU 数据
    auto imu_sub = nh->create_subscription<sensor_msgs::msg::Imu>(
        "/IMU_data", qos, imu_callback);
    
    // 订阅欧拉角数据（如果 talker 中 euler_switch 为 true）
    auto euler_sub = nh->create_subscription<geometry_msgs::msg::Vector3Stamped>(
        "/euler_data", qos, euler_callback);
    
    // 订阅磁力计数据（如果 talker 中 magnetic_switch 为 true）
    auto magnetic_sub = nh->create_subscription<sensor_msgs::msg::MagneticField>(
        "/magnetic_data", qos, magnetic_callback);

    RCLCPP_INFO(nh->get_logger(), "Subscriber started. Waiting for data on /IMU_data, /euler_data, /magnetic_data...");

    rclcpp::spin(nh);
    rclcpp::shutdown();

    return 0;
}
