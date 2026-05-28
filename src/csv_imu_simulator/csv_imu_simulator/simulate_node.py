#!/usr/bin/env python3
import rclpy
from rclpy.node import Node
from sensor_msgs.msg import Imu, MagneticField
from geometry_msgs.msg import Vector3Stamped
import pandas as pd
import time
import math

# 单位转换常数
DEG_TO_RAD = 0.01745329
GRA_ACC = 9.8
UTESLA_TO_TESLA = 0.000001

class CSVImuSimulator(Node):
    def __init__(self):
        super().__init__('csv_imu_simulator')
        
        # 参数
        self.declare_parameter('csv_file', '/home/zhang/products-master/examples/ROS2/data/example_data.csv')
        self.declare_parameter('frame_id', 'base_link')
        self.declare_parameter('imu_topic', '/IMU_data')
        self.declare_parameter('euler_topic', '/euler_data')
        self.declare_parameter('mag_topic', '/magnetic_data')
        self.declare_parameter('publish_imu', True)
        self.declare_parameter('publish_euler', True)
        self.declare_parameter('publish_mag', True)
        
        csv_file = self.get_parameter('csv_file').value
        self.frame_id = self.get_parameter('frame_id').value
        imu_topic = self.get_parameter('imu_topic').value
        euler_topic = self.get_parameter('euler_topic').value
        mag_topic = self.get_parameter('mag_topic').value
        self.pub_imu = self.get_parameter('publish_imu').value
        self.pub_euler = self.get_parameter('publish_euler').value
        self.pub_mag = self.get_parameter('publish_mag').value
        
        # 发布者
        if self.pub_imu:
            self.imu_pub = self.create_publisher(Imu, imu_topic, 10)
        if self.pub_euler:
            self.euler_pub = self.create_publisher(Vector3Stamped, euler_topic, 10)
        if self.pub_mag:
            self.mag_pub = self.create_publisher(MagneticField, mag_topic, 10)
        
        # 读取 CSV，只保留 HI91 行（MATLAB 代码中就是筛选 HI91）
        self.df = pd.read_csv(csv_file)
        # 如果第一列是 frame_type，则筛选；否则假设文件已经只有 HI91 数据
        if 'frame_type' in self.df.columns:
            self.df = self.df[self.df['frame_type'] == 'HI91']
        self.df = self.df.reset_index(drop=True)
        if len(self.df) == 0:
            self.get_logger().error('No HI91 data found in CSV!')
            raise RuntimeError('No HI91 data')
        
        # 确保时间列存在并转换为秒
        if 'sys_time' not in self.df.columns:
            self.get_logger().error('CSV missing sys_time column')
            raise RuntimeError('Missing sys_time')
        # sys_time 可能是整数毫秒，转换为相对秒数
        start_time = self.df['sys_time'].iloc[0] / 1000.0
        self.timestamps_sec = self.df['sys_time'] / 1000.0 - start_time
        
        # 计算帧间隔（用于 sleep）
        self.delays = []
        for i in range(1, len(self.timestamps_sec)):
            dt = self.timestamps_sec.iloc[i] - self.timestamps_sec.iloc[i-1]
            if dt > 0:
                self.delays.append(dt)
        if not self.delays:
            self.delays = [0.01]  # fallback 100Hz
        
        self.idx = 0
        self.start_real_time = None
        self.timer = self.create_timer(0.001, self.publish_next)  # 高频检查，不准确，改用线程循环？更好用定时器重调度
        
        # 改用 while 循环在单独线程，或者使用定时器且每个周期检查时间
        # 更简单的做法：在回调中基于真实时间判断
        self.get_logger().info(f'CSV IMU Simulator started, {len(self.df)} frames')
        
    def publish_next(self):
        if self.idx >= len(self.df):
            self.get_logger().info('Finished publishing all frames.')
            self.destroy_timer(self.timer)
            return
        
        # 检查是否到达发布时间
        if self.start_real_time is None:
            self.start_real_time = self.get_clock().now().nanoseconds / 1e9
        
        target_rel_time = self.timestamps_sec.iloc[self.idx]
        now_rel = (self.get_clock().now().nanoseconds / 1e9) - self.start_real_time
        
        if now_rel >= target_rel_time:
            # 发布当前帧
            row = self.df.iloc[self.idx]
            self.publish_frame(row)
            self.idx += 1
            # 如果已发布完，停止定时器
            if self.idx >= len(self.df):
                self.destroy_timer(self.timer)
                
    def publish_frame(self, row):
        stamp = self.get_clock().now().to_msg()
        
        # 1) IMU 消息
        if self.pub_imu:
            imu_msg = Imu()
            imu_msg.header.stamp = stamp
            imu_msg.header.frame_id = self.frame_id
            
            # 四元数
            imu_msg.orientation.w = float(row['quat_w'])
            imu_msg.orientation.x = float(row['quat_x'])
            imu_msg.orientation.y = float(row['quat_y'])
            imu_msg.orientation.z = float(row['quat_z'])
            # 角速度 (deg/s -> rad/s)
            imu_msg.angular_velocity.x = float(row['gyr_x']) * DEG_TO_RAD
            imu_msg.angular_velocity.y = float(row['gyr_y']) * DEG_TO_RAD
            imu_msg.angular_velocity.z = float(row['gyr_z']) * DEG_TO_RAD
            # 线加速度 (g -> m/s^2)
            imu_msg.linear_acceleration.x = float(row['acc_x']) * GRA_ACC
            imu_msg.linear_acceleration.y = float(row['acc_y']) * GRA_ACC
            imu_msg.linear_acceleration.z = float(row['acc_z']) * GRA_ACC
            
            # 可选：设置协方差（这里用零）
            imu_msg.orientation_covariance[0] = -1.0  # 表示没有协方差
            imu_msg.angular_velocity_covariance[0] = -1.0
            imu_msg.linear_acceleration_covariance[0] = -1.0
            
            self.imu_pub.publish(imu_msg)
        
        # 2) 欧拉角消息
        if self.pub_euler:
            euler_msg = Vector3Stamped()
            euler_msg.header.stamp = stamp
            euler_msg.header.frame_id = self.frame_id
            euler_msg.vector.x = float(row['roll']) * DEG_TO_RAD
            euler_msg.vector.y = float(row['pitch']) * DEG_TO_RAD
            euler_msg.vector.z = float(row['yaw']) * DEG_TO_RAD
            self.euler_pub.publish(euler_msg)
        
        # 3) 磁场消息
        if self.pub_mag:
            mag_msg = MagneticField()
            mag_msg.header.stamp = stamp
            mag_msg.header.frame_id = self.frame_id
            mag_msg.magnetic_field.x = float(row['mag_x']) * UTESLA_TO_TESLA
            mag_msg.magnetic_field.y = float(row['mag_y']) * UTESLA_TO_TESLA
            mag_msg.magnetic_field.z = float(row['mag_z']) * UTESLA_TO_TESLA
            self.mag_pub.publish(mag_msg)
        
        self.get_logger().debug(f'Published frame {self.idx}/{len(self.df)}', throttle_duration_sec=1.0)

def main(args=None):
    rclpy.init(args=args)
    node = CSVImuSimulator()
    rclpy.spin(node)
    node.destroy_node()
    rclpy.shutdown()

if __name__ == '__main__':
    main()