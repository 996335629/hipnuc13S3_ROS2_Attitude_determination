cmake_minimum_required(VERSION 3.10)
project(imu_9dof_attitude)

find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(std_srvs REQUIRED)
find_package(tf2 REQUIRED)
find_package(tf2_geometry_msgs REQUIRED)
find_package(tf2_ros REQUIRED)

add_executable(imu_9dof_attitude src/imu_9dof_attitude.cpp)
target_compile_features(imu_9dof_attitude PUBLIC cxx_std_17)
ament_target_dependencies(imu_9dof_attitude
  rclcpp
  sensor_msgs
  geometry_msgs
  std_srvs
  tf2
  tf2_geometry_msgs
  tf2_ros
)

install(TARGETS
  imu_9dof_attitude
  DESTINATION lib/${PROJECT_NAME}
)

# 安装 launch 文件（关键）
install(DIRECTORY launch
  DESTINATION share/${PROJECT_NAME}
)

ament_package()