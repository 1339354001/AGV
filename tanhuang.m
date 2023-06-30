close all; clc; clear;
% 定义系统的常数
g = 9.81; % 重力加速度(m/s^2)
m = 1; % 摆球质量(kg)
k = 100; % 劲度系数(N/m)
L0 = 1; % 摆长(m)
v0 = 2; % 小球在平衡位置处的速度(m/s)

% 计算初始情况
r0 = L0 + m * g / k;
theta0 = 0;
dr0 = 0;
dtheta0 = v0 / r0;

% 定义时间区间和单位区间长度
tspan = [0, 100];
dt = 0.01;

% 输入初始情况
y0 = [r0; theta0; dr0; dtheta0];

% 定义常微分方程
f = @(t, y) [y(3); y(4); -k/m*(y(1)-L0-m*g/k)-y(1)*y(4)^2+g*cos(y(2)); -2*y(3)*y(4)/y(1)-g*sin(y(2))/y(1)];
% y(1)为摆长，y(2)为角度，y(3)为法向速度，y(4)为角速度

% 进行4阶龙格-库塔法数值微分求解
[t, y] = rk4(f, tspan, y0, dt);

% 提取摆长、法向速度、角度、角速度和切向速度
r = y(:, 1);
theta = y(:, 2);
vn = y(:, 3);
w = y(:, 4);
vt = w.*r;

% 计算直角坐标
x = r .* sin(theta);
y = -r .* cos(theta);

% 输出摆长随时间的图像
figure;
subplot(2, 1, 1);
plot(t, r);
xlabel('Time (s)');
ylabel('Pendulum Length (m)');
title('Pendulum Length - Time');

% 输出法向线速度随时间的图像
subplot(2, 1, 2);
plot(t, vn);
xlabel('Time (s)');
ylabel('Normal Velocity (m/s)');
title('Normal Velocity - Time');

% 输出角度随时间的图像
figure;
subplot(2, 1, 1);
plot(t, theta);
xlabel('Time (s)');
ylabel('Angle (rad)');
title('Angle - Time');

% 输出角速度随时间的图像
subplot(2, 1, 2);
plot(t, w);
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
title('Angular Velocity - Time');

% 输出切向速度随时间的图像
figure;
subplot(2, 1, 1);
plot(t, vt);
xlabel('Time (s)');
ylabel('Tangential Velocity (m/s)');
title('Tangential Velocity - Time');

% 输出摆球运动轨迹
figure;
plot(x, y);
xlabel('Horizontal Position (m)');
ylabel('Vertical Position (m)');
title('Trajectory of Spring Pendulum');

% 龙格-库塔法
function [t, y] = rk4(f, tspan, y0, dt)
    t = tspan(1):dt:tspan(2);% 离散化的输出区间
    n = length(t);% 区间内的小区间个数
    y = zeros(n, length(y0));% 创建输出结果y
    y(1, :) = y0;% y的第1行为初始条件

    % 4阶龙格库塔计算
    for i = 1:(n-1)
        % 计算四个参数
        k1 = f(t(i), y(i, :))';
        k2 = f(t(i) + 0.5*dt, y(i, :) + 0.5*dt*k1)';
        k3 = f(t(i) + 0.5*dt, y(i, :) + 0.5*dt*k2)';
        k4 = f(t(i) + dt, y(i, :) + dt*k3)';
        % 更新下一个y
        y(i+1, :) = y(i, :) + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*dt;
    end
end
