%% indi_attitude_control_sim.m
% 基于 INDI 的滑翔飞行器 BTT 姿态控制仿真
%
% 本版本重点保证动力学结构的完整性：
% 1. 平动状态与转动状态分开建模；
% 2. 航迹倾角 gamma、航向角 chi 单独积分，不再用 theta-alpha 代替；
% 3. 姿态角 theta、phi、psi 由角速度积分得到；
% 4. alpha、beta 用相对低阶气动动态描述，但与 gamma_dot、chi_dot 耦合，
%    使其不再脱离平动动力学单独漂移；
% 5. 采用 BTT 转弯控制：beta 目标为 0，侧向过载主要由滚转后的升力分解建立。
%
% 当前版本仍然是“姿态控制验证模型”，因此没有加入：
% - 对称增阻模态
% - 时间约束控制
% - 舵机动态与舵速限制

clear; clc; close all;

%% -------------------- 参数定义 --------------------
P = struct();

% 常量与飞行器参数
P.g               = 9.81;                 % 重力加速度 m/s^2
P.m               = 500;                  % 质量 kg
P.rho             = 0.90;                 % 空气密度 kg/m^3
P.S               = 0.85;                 % 参考面积 m^2

% 转动惯量
P.Jx              = 120;
P.Jy              = 185;
P.Jz              = 205;
P.J               = diag([P.Jx, P.Jy, P.Jz]);

% 升力、阻力、侧力模型
P.CL0             = 0.05;
P.CL_alpha        = 4.6;
P.CD0             = 0.05;
P.K_induced       = 0.09;
P.CY_beta         = -0.55;

% 力矩模型
P.Lp              = 90;
P.Mq              = 110;
P.Nr              = 85;
P.M_alpha         = 420;
P.N_beta          = 280;
P.L_beta          = 20;

% alpha / beta 动态附加项
P.a_alpha         = 1.6;
P.a_beta          = 2.0;
P.b_alpha         = 1.2;
P.b_beta          = 1.2;

% X 布局尾舵模态分解
P.Ta = [ 1,  1, -1;
        -1,  1,  1;
         1, -1,  1;
        -1, -1, -1];

% 四舵控制效能矩阵
P.Kr_fin          = 120;
P.Kp_fin          = 160;
P.Ky_fin          = 130;
P.B4 = [ P.Kr_fin, -P.Kr_fin,  P.Kr_fin, -P.Kr_fin;
          P.Kp_fin,  P.Kp_fin, -P.Kp_fin, -P.Kp_fin;
         -P.Ky_fin,  P.Ky_fin,  P.Ky_fin, -P.Ky_fin];
P.Ba              = P.B4 * P.Ta;
P.Ba_pinv         = pinv(P.Ba);

% BTT 外环参数
P.k_alpha_q       = 7.5;
P.k_phi_p         = 4.8;
P.k_phi_d         = 2.0;
P.k_beta_r        = 7.5;
P.small_load      = 0.2;
P.phi_cmd_max     = deg2rad(55);
P.alpha_cmd_lim   = deg2rad([-2; 12]);
P.beta_cmd_lim    = deg2rad([-5; 5]);
P.p_cmd_max       = deg2rad(90);
P.q_cmd_max       = deg2rad(70);
P.r_cmd_max       = deg2rad(70);

% INDI 内环参数
P.Komega          = diag([12, 18, 16]);
P.tau_cmd_dot     = 0.05;
P.tau_omegadot    = 0.05;

% 限幅
P.ua_max          = deg2rad([18; 20; 20]);
P.delta_max       = deg2rad(25);
P.alpha_lim       = deg2rad([-4; 14]);
P.beta_lim        = deg2rad([-6; 6]);
P.theta_lim       = deg2rad([-25; 25]);
P.phi_lim         = deg2rad([-85; 85]);
P.gamma_lim       = deg2rad([-30; 8]);
P.V_min           = 80;
P.V_max           = 320;

% 仿真设置
P.dt              = 0.01;
P.T_end           = 50;
P.N               = round(P.T_end / P.dt) + 1;
P.t               = (0:P.N-1) * P.dt;

%% -------------------- 初始条件 --------------------
X = struct();
X.x               = zeros(1, P.N);
X.y               = zeros(1, P.N);
X.z               = zeros(1, P.N);
X.V               = zeros(1, P.N);
X.gamma           = zeros(1, P.N);
X.chi             = zeros(1, P.N);
X.phi             = zeros(1, P.N);
X.theta           = zeros(1, P.N);
X.psi             = zeros(1, P.N);
X.p               = zeros(1, P.N);
X.q               = zeros(1, P.N);
X.r               = zeros(1, P.N);
X.alpha           = zeros(1, P.N);
X.beta            = zeros(1, P.N);

X.z(1)            = 8000;
X.V(1)            = 250;
X.gamma(1)        = deg2rad(-5.0);
X.chi(1)          = 0.0;
X.phi(1)          = 0.0;
X.alpha(1)        = deg2rad(2.0);
X.theta(1)        = X.gamma(1) + X.alpha(1);
X.beta(1)         = 0.0;
X.psi(1)          = X.chi(1);

%% -------------------- 结果存储 --------------------
R = struct();
R.n_n_cmd         = zeros(1, P.N);
R.n_l_cmd         = zeros(1, P.N);
R.n_total_cmd     = zeros(1, P.N);
R.n_n             = zeros(1, P.N);
R.n_l             = zeros(1, P.N);
R.n_total         = zeros(1, P.N);
R.alpha_cmd       = zeros(1, P.N);
R.beta_cmd        = zeros(1, P.N);
R.phi_cmd         = zeros(1, P.N);
R.p_cmd           = zeros(1, P.N);
R.q_cmd           = zeros(1, P.N);
R.r_cmd           = zeros(1, P.N);
R.omega_cmd       = zeros(3, P.N);
R.omega_dot_meas  = zeros(3, P.N);
R.nu              = zeros(3, P.N);
R.ua              = zeros(3, P.N);
R.delta           = zeros(4, P.N);
R.M_ctrl          = zeros(3, P.N);
R.M_aero          = zeros(3, P.N);
R.CL              = zeros(1, P.N);
R.CD              = zeros(1, P.N);
R.CY              = zeros(1, P.N);

[R.n_n(1), R.n_l(1), R.n_total(1), R.CL(1), R.CD(1), R.CY(1)] = ...
    computeLoadsAndCoeffs(X.alpha(1), X.beta(1), X.phi(1), X.V(1), P);

%% -------------------- 主仿真循环 --------------------
ua_prev = zeros(3, 1);
omega_prev = [X.p(1); X.q(1); X.r(1)];
omega_dot_f_prev = zeros(3, 1);
omega_cmd_prev = zeros(3, 1);

for k = 1:P.N-1
    t = P.t(k);

    % 当前状态
    xk = [X.x(k); X.y(k); X.z(k); X.V(k); X.gamma(k); X.chi(k); ...
          X.phi(k); X.theta(k); X.psi(k); X.p(k); X.q(k); X.r(k); ...
          X.alpha(k); X.beta(k)];

    V     = X.V(k);
    phi   = X.phi(k);
    alpha = X.alpha(k);
    beta  = X.beta(k);
    omega = [X.p(k); X.q(k); X.r(k)];

    % 当前动压与总过载增益
    qbar = 0.5 * P.rho * V^2;
    nn_alpha_gain = qbar * P.S * P.CL_alpha / (P.m * P.g);

    % 时变过载指令
    [n_n_cmd, n_l_cmd] = commandProfile(t);
    n_total_cmd = sqrt(max(n_n_cmd, 0)^2 + n_l_cmd^2);

    % BTT 外环：
    % alpha 控制总过载大小，phi 控制过载方向，beta 压向 0。
    alpha_cmd = n_total_cmd / max(nn_alpha_gain, 1e-6);
    alpha_cmd = saturate(alpha_cmd, P.alpha_cmd_lim(1), P.alpha_cmd_lim(2));

    phi_cmd = atan2(n_l_cmd, max(n_n_cmd, P.small_load));
    phi_cmd = saturate(phi_cmd, -P.phi_cmd_max, P.phi_cmd_max);

    beta_cmd = 0.0;

    % 角速度指令
    p_cmd = P.k_phi_p * (phi_cmd - X.phi(k)) - P.k_phi_d * X.p(k);
    q_cmd = P.k_alpha_q * (alpha_cmd - alpha);
    r_cmd = P.k_beta_r  * (beta_cmd  - beta);

    p_cmd = saturate(p_cmd, -P.p_cmd_max, P.p_cmd_max);
    q_cmd = saturate(q_cmd, -P.q_cmd_max, P.q_cmd_max);
    r_cmd = saturate(r_cmd, -P.r_cmd_max, P.r_cmd_max);
    omega_cmd = [p_cmd; q_cmd; r_cmd];

    % 期望角加速度
    omega_cmd_dot = (omega_cmd - omega_cmd_prev) / max(P.tau_cmd_dot, P.dt);
    nu = omega_cmd_dot + P.Komega * (omega_cmd - omega);

    % 角加速度估计
    omega_dot_raw = (omega - omega_prev) / P.dt;
    omega_dot_meas = omega_dot_f_prev + ...
        (P.dt / P.tau_omegadot) * (omega_dot_raw - omega_dot_f_prev);

    % INDI 控制律
    Delta_ua = P.Ba_pinv * (P.J * (nu - omega_dot_f_prev));
    ua_cmd = ua_prev + Delta_ua;
    ua_cmd = min(max(ua_cmd, -P.ua_max), P.ua_max);

    % 模态到四舵映射并限幅
    [delta_cmd, ua_eff] = mixAndLimitXtail(ua_cmd, P);
    M_ctrl = P.B4 * delta_cmd;

    uk.ua = ua_eff;
    uk.M_ctrl = M_ctrl;

    % RK4 积分
    xk1 = rk4Step(@(x) plantDynamics(x, uk, P), xk, P.dt);

    % 状态更新
    X.x(k+1)      = xk1(1);
    X.y(k+1)      = xk1(2);
    X.z(k+1)      = max(xk1(3), 0);
    X.V(k+1)      = saturate(xk1(4), P.V_min, P.V_max);
    X.gamma(k+1)  = saturate(xk1(5), P.gamma_lim(1), P.gamma_lim(2));
    X.chi(k+1)    = wrapToPiLocal(xk1(6));
    X.phi(k+1)    = saturate(wrapToPiLocal(xk1(7)), P.phi_lim(1), P.phi_lim(2));
    X.theta(k+1)  = saturate(xk1(8), P.theta_lim(1), P.theta_lim(2));
    X.psi(k+1)    = wrapToPiLocal(xk1(9));
    X.p(k+1)      = xk1(10);
    X.q(k+1)      = xk1(11);
    X.r(k+1)      = xk1(12);
    X.alpha(k+1)  = saturate(xk1(13), P.alpha_lim(1), P.alpha_lim(2));
    X.beta(k+1)   = saturate(xk1(14), P.beta_lim(1), P.beta_lim(2));

    [n_n, n_l, n_total, CL, CD, CY] = ...
        computeLoadsAndCoeffs(X.alpha(k+1), X.beta(k+1), X.phi(k+1), X.V(k+1), P);
    M_aero = aeroMoments([X.p(k); X.q(k); X.r(k)], X.alpha(k), X.beta(k), P);

    % 记录
    R.n_n_cmd(k)          = n_n_cmd;
    R.n_l_cmd(k)          = n_l_cmd;
    R.n_total_cmd(k)      = n_total_cmd;
    R.n_n(k+1)            = n_n;
    R.n_l(k+1)            = n_l;
    R.n_total(k+1)        = n_total;
    R.alpha_cmd(k)        = alpha_cmd;
    R.beta_cmd(k)         = beta_cmd;
    R.phi_cmd(k)          = phi_cmd;
    R.p_cmd(k)            = p_cmd;
    R.q_cmd(k)            = q_cmd;
    R.r_cmd(k)            = r_cmd;
    R.omega_cmd(:, k)     = omega_cmd;
    R.omega_dot_meas(:, k)= omega_dot_meas;
    R.nu(:, k)            = nu;
    R.ua(:, k)            = ua_eff;
    R.delta(:, k)         = delta_cmd;
    R.M_ctrl(:, k)        = M_ctrl;
    R.M_aero(:, k)        = M_aero;
    R.CL(k+1)             = CL;
    R.CD(k+1)             = CD;
    R.CY(k+1)             = CY;

    % 暂存更新
    ua_prev = ua_eff;
    omega_prev = omega;
    omega_dot_f_prev = omega_dot_meas;
    omega_cmd_prev = omega_cmd;
end

% 末端补齐
[R.n_n_cmd(end), R.n_l_cmd(end)] = commandProfile(P.t(end));
R.n_total_cmd(end) = sqrt(max(R.n_n_cmd(end), 0)^2 + R.n_l_cmd(end)^2);
R.alpha_cmd(end) = R.alpha_cmd(end-1);
R.beta_cmd(end)  = R.beta_cmd(end-1);
R.phi_cmd(end)   = R.phi_cmd(end-1);
R.p_cmd(end)     = R.p_cmd(end-1);
R.q_cmd(end)     = R.q_cmd(end-1);
R.r_cmd(end)     = R.r_cmd(end-1);
R.omega_cmd(:, end) = R.omega_cmd(:, end-1);
R.omega_dot_meas(:, end) = R.omega_dot_meas(:, end-1);
R.nu(:, end) = R.nu(:, end-1);
R.ua(:, end) = R.ua(:, end-1);
R.delta(:, end) = R.delta(:, end-1);
R.M_ctrl(:, end) = R.M_ctrl(:, end-1);
R.M_aero(:, end) = R.M_aero(:, end-1);

%% -------------------- 绘图 --------------------
t = P.t;
deg = 180 / pi;

% 1. 弹道曲线
figure(1); clf;
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(X.x / 1000, X.z / 1000, 'b', 'LineWidth', 1.8);
grid on;
xlabel('x (km)');
ylabel('z (km)');
title('纵向弹道曲线');
nexttile;
plot3(X.x / 1000, X.y / 1000, X.z / 1000, 'r', 'LineWidth', 1.8);
grid on;
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
title('三维弹道曲线');
view(35, 24);

% 2. 速度曲线
figure(2); clf;
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t, X.V, 'b', 'LineWidth', 1.8);
grid on;
xlabel('时间 (s)');
ylabel('V (m/s)');
title('速度曲线');
nexttile;
plot(t, X.gamma * deg, 'm', 'LineWidth', 1.8);
grid on;
xlabel('时间 (s)');
ylabel('\gamma (deg)');
title('航迹倾角');

% 3. 姿态曲线
figure(3); clf;
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t, X.phi * deg, 'b', 'LineWidth', 1.8); hold on;
plot(t, R.phi_cmd * deg, 'k--', 'LineWidth', 1.3);
grid on;
xlabel('时间 (s)');
ylabel('\phi (deg)');
legend('\phi', '\phi_{cmd}', 'Location', 'best');
title('滚转角');
nexttile;
plot(t, X.theta * deg, 'b', 'LineWidth', 1.8);
grid on;
xlabel('时间 (s)');
ylabel('\theta (deg)');
title('俯仰角');
nexttile;
plot(t, X.psi * deg, 'b', 'LineWidth', 1.8);
grid on;
xlabel('时间 (s)');
ylabel('\psi (deg)');
title('航向角');

% 4. 过载跟踪曲线
figure(4); clf;
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t, R.n_n_cmd, 'k--', 'LineWidth', 1.4); hold on;
plot(t, R.n_n, 'b', 'LineWidth', 1.8);
grid on;
xlabel('时间 (s)');
ylabel('n_n (g)');
legend('n_{n,cmd}', 'n_n', 'Location', 'best');
title('法向过载跟踪');
nexttile;
plot(t, R.n_l_cmd, 'k--', 'LineWidth', 1.4); hold on;
plot(t, R.n_l, 'r', 'LineWidth', 1.8);
grid on;
xlabel('时间 (s)');
ylabel('n_l (g)');
legend('n_{l,cmd}', 'n_l', 'Location', 'best');
title('侧向过载跟踪');
nexttile;
plot(t, R.n_total_cmd, 'k--', 'LineWidth', 1.4); hold on;
plot(t, R.n_total, 'm', 'LineWidth', 1.8);
grid on;
xlabel('时间 (s)');
ylabel('n_{tot} (g)');
legend('n_{tot,cmd}', 'n_{tot}', 'Location', 'best');
title('总过载');

% 5. 角度跟踪曲线
figure(5); clf;
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t, R.alpha_cmd * deg, 'k--', 'LineWidth', 1.4); hold on;
plot(t, X.alpha * deg, 'b', 'LineWidth', 1.8);
grid on;
xlabel('时间 (s)');
ylabel('\alpha (deg)');
legend('\alpha_{cmd}', '\alpha', 'Location', 'best');
title('攻角跟踪');
nexttile;
plot(t, R.beta_cmd * deg, 'k--', 'LineWidth', 1.4); hold on;
plot(t, X.beta * deg, 'r', 'LineWidth', 1.8);
grid on;
xlabel('时间 (s)');
ylabel('\beta (deg)');
legend('\beta_{cmd}', '\beta', 'Location', 'best');
title('侧滑角跟踪');
nexttile;
plot(t, R.phi_cmd * deg, 'k--', 'LineWidth', 1.4); hold on;
plot(t, X.phi * deg, 'b', 'LineWidth', 1.8);
grid on;
xlabel('时间 (s)');
ylabel('\phi (deg)');
legend('\phi_{cmd}', '\phi', 'Location', 'best');
title('滚转角跟踪');

% 6. 舵偏角曲线
figure(6); clf;
tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:4
    nexttile;
    plot(t, R.delta(i, :) * deg, 'LineWidth', 1.8);
    grid on;
    xlabel('时间 (s)');
    ylabel(sprintf('\\delta_%d (deg)', i));
    title(sprintf('尾舵偏角 \\delta_%d', i));
end

% 7. 角速度曲线
figure(7); clf;
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t, X.p * deg, 'b', 'LineWidth', 1.8); hold on;
plot(t, R.p_cmd * deg, 'k--', 'LineWidth', 1.3);
grid on;
xlabel('时间 (s)');
ylabel('p (deg/s)');
legend('p', 'p_{cmd}', 'Location', 'best');
title('滚转角速度');
nexttile;
plot(t, X.q * deg, 'b', 'LineWidth', 1.8); hold on;
plot(t, R.q_cmd * deg, 'k--', 'LineWidth', 1.3);
grid on;
xlabel('时间 (s)');
ylabel('q (deg/s)');
legend('q', 'q_{cmd}', 'Location', 'best');
title('俯仰角速度');
nexttile;
plot(t, X.r * deg, 'b', 'LineWidth', 1.8); hold on;
plot(t, R.r_cmd * deg, 'k--', 'LineWidth', 1.3);
grid on;
xlabel('时间 (s)');
ylabel('r (deg/s)');
legend('r', 'r_{cmd}', 'Location', 'best');
title('偏航角速度');

%% -------------------- 结果摘要 --------------------
R.summary = buildSummary(R, X, P);

fprintf('\nBTT-INDI 姿态控制仿真完成。\n');
fprintf('法向过载峰值误差: %.3f g\n', R.summary.peak_nn_error_g);
fprintf('侧向过载峰值误差: %.3f g\n', R.summary.peak_nl_error_g);
fprintf('最大舵偏角: %.2f deg\n', R.summary.max_delta_deg);
fprintf('末端速度: %.2f m/s\n', X.V(end));
fprintf('末端俯仰角: %.2f deg\n', X.theta(end) * deg);
fprintf('末端航迹角: %.2f deg\n', X.gamma(end) * deg);
fprintf('Ba 条件数: %.2f\n', R.summary.cond_Ba);

disp(' ');
disp('当前版本的主要简化：');
disp('- 使用 BTT 转弯控制，beta 目标固定为 0，侧向过载主要通过滚转角建立。');
disp('- 已单独积分 V、gamma、chi、phi、theta、psi、p、q、r、alpha、beta，动力学结构更完整。');
disp('- alpha/beta 仍采用低阶近似模型，因此仍属于控制律验证模型，而非高保真气动模型。');
disp('- 不包含舵机动态、对称增阻模态、时间控制通道。');

%% ==================== 本地函数 ====================
function [n_n_cmd, n_l_cmd] = commandProfile(t)
% 时变过载指令

    if t < 5
        n_n_cmd = 1.0;
        n_l_cmd = 0.0;
    elseif t < 15
        n_n_cmd = 1.0 + 0.12 * (t - 5);
        n_l_cmd = 0.0;
    elseif t < 25
        n_n_cmd = 2.2;
        n_l_cmd = 0.05 * (t - 15);
    elseif t < 35
        n_n_cmd = 2.2 - 0.08 * (t - 25);
        n_l_cmd = 0.5 - 0.13 * (t - 25);
    else
        n_n_cmd = 1.2 + 0.10 * sin(0.65 * (t - 35));
        n_l_cmd = -0.10 + 0.18 * sin(0.45 * (t - 35));
    end
end

function dx = plantDynamics(x, u, P)
% 简化但结构完整的平动+转动动力学
% 状态：
% x = [x; y; z; V; gamma; chi; phi; theta; psi; p; q; r; alpha; beta]

    V     = max(x(4), P.V_min);
    gamma = x(5);
    chi   = x(6);
    phi   = x(7);
    theta = x(8);
    psi   = x(9);
    p     = x(10);
    q     = x(11);
    r     = x(12);
    alpha = x(13);
    beta  = x(14);

    omega = [p; q; r];

    [n_n, n_l, n_total, CL, CD, CY] = computeLoadsAndCoeffs(alpha, beta, phi, V, P); %#ok<ASGLU>
    qbar = 0.5 * P.rho * V^2;
    L = qbar * P.S * CL;
    D = qbar * P.S * CD;
    Y = qbar * P.S * CY;

    % 平动动力学：速度、航迹倾角、航向角分开积分
    x_dot = V * cos(gamma) * cos(chi);
    y_dot = V * cos(gamma) * sin(chi);
    z_dot = V * sin(gamma);
    V_dot = -D / P.m - P.g * sin(gamma);
    gamma_dot = (L * cos(phi) / P.m - P.g * cos(gamma)) / max(V, P.V_min);
    chi_dot = (L * sin(phi) / (P.m * max(V, P.V_min) * max(cos(gamma), 0.08)));

    % 欧拉角运动学
    T = eulerRateMatrix(phi, theta);
    euler_dot = T * omega;

    % 转动动力学
    M_aero = aeroMoments(omega, alpha, beta, P);
    omega_dot = P.J \ (M_aero + u.M_ctrl - cross(omega, P.J * omega));

    % alpha / beta 动态与平动耦合
    % alpha 反映机体俯仰相对航迹倾角的变化，因此减去 gamma_dot
    % beta 反映机体航向相对速度方向的偏差，因此减去 chi_dot
    alpha_dot = q - gamma_dot - P.a_alpha * alpha + P.b_alpha * u.ua(2);
    beta_dot  = r - chi_dot   - P.a_beta  * beta  + P.b_beta  * u.ua(3);

    dx = [x_dot; y_dot; z_dot; V_dot; gamma_dot; chi_dot; ...
          euler_dot; omega_dot; alpha_dot; beta_dot];
end

function M_aero = aeroMoments(omega, alpha, beta, P)
% 简化气动力矩模型

    p = omega(1);
    q = omega(2);
    r = omega(3);

    M_aero = [ -P.Lp * p - P.L_beta * beta;
               -P.Mq * q - P.M_alpha * alpha;
               -P.Nr * r - P.N_beta * beta ];
end

function [n_n, n_l, n_total, CL, CD, CY] = computeLoadsAndCoeffs(alpha, beta, phi, V, P)
% BTT 过载与气动系数模型

    qbar = 0.5 * P.rho * V^2;
    CL = P.CL0 + P.CL_alpha * alpha;
    CD = P.CD0 + P.K_induced * CL^2 + 0.02 * beta^2;
    CY = P.CY_beta * beta;

    L = qbar * P.S * CL;
    Y = qbar * P.S * CY;

    n_total = L / (P.m * P.g);
    n_n = n_total * cos(phi);
    n_l = n_total * sin(phi) + Y / (P.m * P.g);
end

function [delta_limited, ua_eff] = mixAndLimitXtail(ua_cmd, P)
% X 尾翼模态到四舵映射，并做统一缩放限幅

    delta_cmd = P.Ta * ua_cmd;
    max_abs = max(abs(delta_cmd));

    if max_abs > P.delta_max
        scale = P.delta_max / max_abs;
    else
        scale = 1.0;
    end

    delta_limited = delta_cmd * scale;
    ua_eff = pinv(P.Ta) * delta_limited;
    ua_eff = min(max(ua_eff, -P.ua_max), P.ua_max);
end

function T = eulerRateMatrix(phi, theta)
% 欧拉角运动学矩阵

    ct = cos(theta);
    if abs(ct) < 0.08
        ct = sign(ct + eps) * 0.08;
    end

    T = [1, sin(phi) * tan(theta),  cos(phi) * tan(theta);
         0, cos(phi),              -sin(phi);
         0, sin(phi) / ct,          cos(phi) / ct];
end

function x_next = rk4Step(f, x, dt)
% 四阶 Runge-Kutta 单步积分

    k1 = f(x);
    k2 = f(x + 0.5 * dt * k1);
    k3 = f(x + 0.5 * dt * k2);
    k4 = f(x + dt * k3);
    x_next = x + dt * (k1 + 2*k2 + 2*k3 + k4) / 6;
end

function y = saturate(u, umin, umax)
% 限幅函数
    y = min(max(u, umin), umax);
end

function ang = wrapToPiLocal(ang)
% 角度归一化到 [-pi, pi]
    ang = mod(ang + pi, 2*pi) - pi;
end

function summary = buildSummary(R, X, P)
% 结果摘要

    summary = struct();
    summary.peak_nn_error_g = max(abs(R.n_n_cmd - R.n_n));
    summary.peak_nl_error_g = max(abs(R.n_l_cmd - R.n_l));
    summary.rms_nn_error_g  = sqrt(mean((R.n_n_cmd - R.n_n).^2));
    summary.rms_nl_error_g  = sqrt(mean((R.n_l_cmd - R.n_l).^2));
    summary.max_delta_deg   = max(abs(R.delta(:))) * 180 / pi;
    summary.max_modal_deg   = max(abs(R.ua(:))) * 180 / pi;
    summary.final_speed     = X.V(end);
    summary.final_range     = sqrt(X.x(end)^2 + X.y(end)^2);
    summary.cond_Ba         = cond(P.Ba);
end
