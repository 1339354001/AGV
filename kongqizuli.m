close all; clc; clear;
% 定义系统的常数
g = 9.81; % 重力加速度(m/s^2)
l = 1; % 摆线长度(m)
m = 1; % 摆球质量(kg)
R = 0.1; % 摆球半径(m)
rho = 1.2; % 空气密度(kg/m^3)
C = 2; % 空气阻力系数
theta0 = pi/4; % 初始摆角(rad)
omega0 = 0; % 初始角速度(rad/s)

% 定义时间区间和单位区间长度
tspan = [0 100];
dt = 0.01;

% 输入初始情况
y0 = [theta0; omega0];

% 定义常微分方程   
f = @(t, y) [y(2); -g/l*sin(y(1)) - 0.5*rho*C*pi*R^2/m*norm(y(2))*y(2)*(l+R)];% y(1)为角度，y(2)为角速度

% 进行4阶龙格-库塔法数值微分求解
[t, y] = rk4(f, tspan, y0, dt);

% 提取角度和角速度
theta = y(:, 1);
omega = y(:, 2);

% 计算直角坐标
x = l*sin(theta);
y = -l*cos(theta);

% 可视化
% 输出角度随时间的图像
subplot(2, 1, 1);
plot(t, theta);
xlabel('Time (s)');
ylabel('Angle (rad)');
title('Angle - Time');

% 输出角速度随时间的图像
subplot(2, 1, 2);
plot(t, omega);
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
title('Angular Velocity - Time');

% 输出摆球运动轨迹
figure;
plot(x, y);
xlabel('Horizontal Position (m)');
ylabel('Vertical Position (m)');
title('Trajectory of Pendulum Ball');

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