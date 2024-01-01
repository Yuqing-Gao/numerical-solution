clear;clc

% 对空间离散化
xmin = 0;
xmax = 1;
ymin = 0;
ymax = 1;
N = 20;
dx = (xmax - xmin) / (N - 1); % 步长
dy = (ymax - ymin) / (N - 1);
[x, y] = meshgrid(xmin:dx:xmax, ymin:dy:ymax);

% 对时间离散化
tmin = 0;
tmax = 1;
t_steps = 100;
dt = (tmax - tmin) / t_steps; % 时间步长
t = tmin:dt:tmax;

% 计算隐式方程的矩阵
A = zeros(N^2, N^2);
B = zeros(N^2, N^2);

% 系数mu
mu = dt / (dx^2);

% 初值
u(:, :, 1) = sin(pi*x).*sin(pi*y);

% 边界条件
u(1, :, :) = 0;
u(:, 1, :) = 0;
u(N+1, :, :) = 0;
u(:, N+1, :) = 0;

for i = 1:N^2
    if i <= N || i >= N^2 - N
        A(i, i) = 1;
        B(i, i) = 1;
    elseif mod(i, N) == 0 || mod(i - 1, N) == 0
        A(i, i) = 1;
        B(i, i) = 1;
    elseif mod(i - 2, N) == 0
        A(i, i) = 1 - 2 * mu;
        A(i, i + 1) = mu / 2;
        B(i, i) = 1 + 2 * mu;
        B(i, i + 1) = -mu / 2;
        A(i, i + N) = mu / 2;
        B(i, i + N) = -mu / 2;
        A(i + N, i) = mu / 2;
        B(i + N, i) = -mu / 2;
    elseif mod(i + 1, N) == 0
        A(i, i) = 1 - 2 * mu;
        A(i, i - 1) = mu / 2;
        B(i, i) = 1 + 2 * mu;
        B(i, i - 1) = -mu / 2;
        A(i, i + N) = mu / 2;
        B(i, i + N) = -mu / 2;
        A(i + N, i) = mu / 2;
        B(i + N, i) = -mu / 2;
    else
        A(i, i) = 1 - 2 * mu;
        A(i, i + 1) = mu / 2;
        A(i, i - 1) = mu / 2;
        B(i, i) = 1 + 2 * mu;
        B(i, i + 1) = -mu / 2;
        B(i, i - 1) = -mu / 2;
        A(i, i + N) = mu / 2;
        B(i, i + N) = -mu / 2;
        A(i + N, i) = mu / 2;
        B(i + N, i) = -mu / 2;
    end
end

for i = N^2 - 2 * N:N^2 - N
    A(i, i + N) = 0;
    B(i, i + N) = 0;
    A(i + N, i) = 0;
    B(i + N, i) = 0;
end

% 储存初值的矩阵
XY = zeros(N);
for j = 1:N
    for i = 1:N
        XY(j, i) = u(i, j, 1);
    end
end

step1 = reshape(XY,[],1);

for n = 2:t_steps
    step2 = B \ A * step1; % 求解 B u_n+1 = A u_n
    step1 = step2;
    
    subplot(1,3,1)
    mesh(reshape(step2, [N N]))
    title(['数值解 t = ',num2str(n/t_steps)])

    u_exact = sin(pi*x).*sin(pi*y).*exp(-2*pi^2*t(n-1));
    subplot(1,3,2)
    mesh(u_exact)
    title("精确解")

    subplot(1,3,3)
    mesh(reshape(step2, [N N])-u_exact)
    title("误差")
    hold off
    pause(0.01)
end
