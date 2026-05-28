# 手动实现IMU姿态解算

1.在连接后IMU的串口后，给与权限，sudo chmod 777 /dev/ttyUSB0
（可以永久设置权限，指令为 sudo usermod -aG dialout username
其中username为用户名，把此用户名加入dialout用户组（dialout是一个group，它主要负责对于串口的权限）

2.在~/products-master/examples/ROS2/hipnuc_ws下运行colcon build构建程序

3.可以在config下的
hipnuc_config.yaml 设置发布的原始数据
imu_6dof_config.yaml 设置算法的各项参数

4.打开第一个集成终端（命令行窗口），install source/setup.bash加载加载当前工作空间的环境变量和路径，
运行ros2 launch hipnuc_imu imu_spec_msg.launch.py 后会通过串口处理后获取原始数据，发布话题（数据），包括角速度，加速度,四元数（/IMU_data），磁力计（/magnetic_data）IMU内置的动态卡尔曼滤波解算算法得到的欧拉角（/euler_data）。

打开第二个集成终端（命令行窗口），install source/setup.bash加载加载当前工作空间的环境变量和路径，，运行ros2 run hipnuc_imu EKF_6dof --ros-args --params-file ./src/hipnuc_imu/config/imu_6dof_config.yaml可以获得通过扩展卡尔曼滤波得到的欧拉角

打开第三个集成终端（命令行窗口），install source/setup.bash加载加载当前工作空间的环境变量和路径，
运行ros2 run hipnuc_imu attitude_estimator_6dof --ros-args --params-file ./src/hipnuc_imu/config/imu_6dof_config.yaml可以获得通过Mahony或Madgwick滤波得到的欧拉角（通过设置config/imu_6dof_config.yaml）

# 其他备注：
假如要移植的话请先参考其他案例下的串口获取数据示例(https://github.com/hipnuc/products.git)

6轴无法保证yaw角的稳定，一定会漂。频率越高更新越快飘得越厉害。

6轴算是完成了，EKF效果还可以，但是9轴的yaw角一直搞不定，初始化都无法收敛到正确角度，甚至是官方的给的数据example_data.csv也对不上，建议使用官方输出的四元数,欧拉角吧。

# 目标
IMU9轴姿态解算

更新频率1000Hz

俯仰、横滚参数

量程		X:±180°  Y:±90°
倾角精度		0.1°
分辨率	水平放置	0.01°

航向角参数

量程		Z:±180°
航向精度		9轴算法，磁场校准，动/静态	1°（不受磁场干扰情况下）【1】
		6轴算法，静态		0.5°（动态存在积分累计误差）【2】
分辨率		水平放置			0.01°

算法：mahony，madgwick，EKF

# IMU的参数
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
    导航到机体旋转矩阵:%与上面转置后是一样的，可以带入数值计算看看
    C_n2B =
[cos(roll)*cos(yaw) - sin(pitch)*sin(roll)*sin(yaw),  cos(roll)*sin(yaw) + cos(yaw)*sin(pitch)*sin(roll), -cos(pitch)*sin(roll)]
[-cos(pitch)*sin(yaw),                                             cos(pitch)*cos(yaw),                                              sin(pitch)]
[cos(yaw)*sin(roll) + cos(roll)*sin(pitch)*sin(yaw), sin(roll)*sin(yaw) - cos(roll)*cos(yaw)*sin(pitch) ,  cos(pitch)*cos(roll)]
    *======================================================================*/

# MATLAB原理推导代码
syms pitch roll yaw H V real
%%定义基本旋转矩阵（坐标变换形式）
Rx = @(a) [1 0 0; 0 cos(a) sin(a); 0 -sin(a) cos(a)];
Ry = @(a) [cos(a) 0 -sin(a); 0 1 0; sin(a) 0 cos(a)];
Rz = @(a) [cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 1];

%%从导航系到机体系
C_n2B = Ry(roll) * Rx(pitch) * Rz(yaw);
%C_n2B的结果与文档的四元数方向余弦矩阵是一样的

%%地磁矢量（导航系）
m_n = [0; H; V];
% 机体系磁场
m_b = simplify(C_n2B * m_n);
% 倾斜补偿（从机体系到水平系）
R_horiz = Rx(-1*pitch) * Ry(-1*roll);%转回水平，获得X和Y的夹角

%syms mbx mby mbz
%m_b=[mbx mby mbz].';%%IMU输出的值就是m_b 测试用 

m_h = simplify(R_horiz * m_b)%
%%m_h=[H*sin(yaw) H*cos(yaw) V].';

%%roll=-atan2(ax, az);
%%pitch = atan2(ay, sqrt(ax * ax + az * az));
%%yaw=arctan(mhx,mhy);



# ROS2串口例程

本文档介绍如何在ROS2下来读取超核电子IMU&GNSS的数据，并提供了c++语言例程代码，通过执行ROS2命令，运行相应的节点，就可以看到打印到终端上的信息。

* 测试环境：Ubuntu20.04   

* ROS版本：ROS2 Foxy

* 测试设备：超核电子IMU系列产品

## 安装USB-UART驱动

Ubuntu 系统自带CP210x的驱动，默认不需要安装串口驱动。将调试版连接到电脑上时，会自动识别设备。识别成功后，会在dev目录下出现一个对应的设备:ttyUSBx

检查USB-UART设备是否被Ubantu识别：

1. 打开终端，输入`ls /dev`,先查看已经存在的串口设备。
2. 查看是否已经存在  ttyUSBx 这个设备文件，便于确认对应的端口号。
4. 接下来插入USB线，连接调试板，然后再次执行`ls /dev`。 dev目录下多了一个设备`ttyUSB0`：

```shell
linux@ubuntu:~$ ls /dev
.....
hpet             net           tty11     tty4   ttyS0      ttyUSB0    vhost-vsock
hugepages        null          tty12     tty40  ttyS1      udmabuf  vmci
......
```

4.打开USB设备的可执行权限：

```shell
   $ sudo chmod 777 /dev/ttyUSB0
```

## 编译hipnuc_ws工作空间

1. 打开终端进入/examples/ROS2/hipnuc_ws 目录
2. 执行`colcon build`命令，编译成功后出现如下信息。

```shell
linux@ubuntu20:~/hipnuc_ws$ colcon build
Starting >>> hipnuc_gnss
Starting >>> hipnuc_imu
Finished <<< hipnuc_imu [0.44s]                                      
Finished <<< hipnuc_gnss [0.49s]

Summary: 2 packages finished [0.61s]
linux@ubuntu20:~/hipnuc_ws$ 
```

##  修改串口波特率和设备号

1. 在Ubuntu环境中，支持的波特率为115200, 460800, 921600。本例程使用的默认波特率是115200，默认打开的串口名称是/dev/ttyUSB0。	

2. 如果您需要更高的输出频率，请修改`config/hipnuc_config.yaml`文件中的配置参数。	

```c
#hipnuc_imu config file
IMU_publisher:
    ros__parameters:
        serial_port: "/dev/ttyUSB0"
        baud_rate: 115200
        frame_id: "base_link"
        imu_topic: "/IMU_data"

            
 #hipnuc_gnss config file
 INS_publisher:
    ros__parameters:
        serial_port: "/dev/ttyUSB0"
        baud_rate: 115200
        frame_id: "gnss_link"
        imu_topic: "/rawimu_data"
        nav_topic: "/NavSatFix_data"
```

注意修改后需要回到hipnuc_ws目录下，重新执行`colcon build`命令

## 显示数据
​	查看数据方式：

​	1、输出ROS 定义的sensor_msgs::Imu。

​	2、输出ROS 定义的sensor_msgs::NavSatFix

###  输出ROS标准 Imu.msg

1. 打开终端，执行：

```shell
linux@ubuntu20:~$ ros2 launch hipnuc_imu imu_spec_msg.launch.py
```

​	2.如果执行失败，提示找不到相应的launch文件，则需要配置环境，在当前终端执行：

```shell
linux@ubuntu:~$source <hipnuc_ws_dir>/install/setup.bash
```

​	3.执行成功后，就可以看到所有的信息：

```c
[listener-2] ---
[listener-2] header:
[listener-2] 	stamp:
[listener-2] 	  secs:1639099575
[listener-2] 	  nanosecs:538349240
[listener-2] 	frame_id:base_link
[listener-2] orientation:
[listener-2] 	x: -0.095125280320644379
[listener-2] 	y: -0.483648955821990967
[listener-2] 	z: 0.053129896521568298
[listener-2] 	w: 0.868453860282897949
[listener-2] orientation_covariance: [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
[listener-2] angular_velocity: 
[listener-2] 	x: -0.000815955184543841
[listener-2] 	y: -0.001057390143056437
[listener-2] 	z: 0.001062464062371403
[listener-2] angular_velocity_covariance: [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
[listener-2] linear_acceleration:
[listener-2] 	x: 8.110355603694916482
[listener-2] 	y: -2.125157430768013000
[listener-2] 	z: 5.013053989410400924
[listener-2] linear_acceleration_covariance: [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
[listener-2] ---
```

​	4、另开一个终端窗口，执行`ros2 topic hz /Imu_data`，可以查看话题发布的频率。

```shell
linux@ubuntu20:~$ ros2 topic hz /Imu_data 
average rate: 100.032
	min: 0.008s max: 0.012s std dev: 0.00058s window: 102
average rate: 100.014
	min: 0.008s max: 0.012s std dev: 0.00054s window: 202
average rate: 100.019
	min: 0.007s max: 0.013s std dev: 0.00064s window: 303
^C
linux@ubuntu20:~$ 
```

### 输出ROS标准的NavSatFix.msg

​	1、打开终端，执行：

```shell
linux@ubuntu20:~$ ros2 launch hipnuc_gnss nav_spec_msg.launch.py 
```

​	2、如果执行失败，提示找不到相应的launch文件，则需要配置环境，在当前终端执行：

```shell
linux@ubuntu:~$source <hipnuc_ws_dir>/install/setup.bash
```

​	3、执行成功后，可以看到如下信息：

```shell
[listener_INS-2] header:
[listener_INS-2] 	stamp:
[listener_INS-2] 	  secs: 1724034005
[listener_INS-2] 	  nanosecs: 370900173
[listener_INS-2] 	frame_id: gnss_link
[listener_INS-2] status:
[listener_INS-2] 	status: 1
[listener_INS-2] 	service: 0
[listener_INS-2] latitude: 40.20336080
[listener_INS-2] longitude: 116.24086010
[listener_INS-2] altitude: 66.30100000
[listener_INS-2] orientation_covariance: [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
[listener_INS-2] position_covariance_type: 0
```

​	4、另开一个终端窗口，执行`ros2 topic hz /NavSatFix_data`，可以查看话题发布的频率。

```shell
linux@ubuntu20:~$ ros2 topic hz /NavSatFix_data`
average rate: 10.032
	min: 0.008s max: 0.012s std dev: 0.00058s window: 10
average rate: 10.014
	min: 0.008s max: 0.012s std dev: 0.00054s window: 20
average rate: 10.019
	min: 0.007s max: 0.013s std dev: 0.00064s window: 30
^C
linux@ubuntu20:~$
```

