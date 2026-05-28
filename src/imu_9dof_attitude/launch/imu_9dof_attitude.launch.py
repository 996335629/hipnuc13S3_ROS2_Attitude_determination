from launch import LaunchDescription
from launch_ros.actions import Node

def generate_launch_description():
    return LaunchDescription([
        Node(
            package='imu_9dof_attitude',
            executable='imu_9dof_attitude',
            name='imu_attitude',
            parameters=[{
                'use_mag': False,        # 先禁用磁力计（避免磁场干扰）
                'algorithm': 'ekf',   # 可选 mahony, madgwick, ekf
                'sample_rate': 1000.0,
            }],
            output='screen',
            emulate_tty=True,
        )
    ])