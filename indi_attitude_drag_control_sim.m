%% indi_attitude_drag_control_sim.m
% 基于 INDI 的滑翔飞行器 BTT 姿态 + 阻力联合控制仿真
%
% 本版本的主要改动：
% 1. 在原有 X 尾翼姿态模态基础上增加对称增阻模态 u_d。
% 2. 将原来的时变过载测试改为恒定法向过载指令和恒定速度指令。
% 3. 保持姿态控制为主任务，阻力控制仅使用尾舵剩余控制裕度。

clear; clc; close all;

%% -------------------- 参数定义 --------------------
P = struct();

% 飞行器与环境参数
P.g               = 9.81;                 % m/s^2
P.m               = 500;                  % kg
P.rho             = 0.90;                 % kg/m^3
P.S               = 0.85;                 % m^2

% 转动惯量
P.Jx              = 120;
P.Jy              = 185;
P.Jz              = 205;
P.J               = diag([P.Jx, P.Jy, P.Jz]);

% 气动参数
P.CL0             = 0.05;
P.CL_alpha        = 4.6;
P.CD0             = 0.05;
P.K_induced       = 0.09;
P.CD_beta         = 0.02;
P.CD_ud           = 0.32;                 % 对称舵偏引入的附加阻力
P.CY_beta         = -0.55;

% 气动力矩参数
P.Lp              = 90;
P.Mq              = 110;
P.Nr              = 85;
P.M_alpha         = 420;
P.N_beta          = 280;
P.L_beta          = 20;

% alpha / beta 低阶动态参数
P.a_alpha         = 1.6;
P.a_beta          = 2.0;
P.b_alpha         = 1.2;
P.b_beta          = 1.2;

% X 尾翼模态分解
P.Ta = [ 1,  1, -1;
        -1,  1,  1;
         1, -1,  1;
        -1, -1, -1];
P.Td              = ones(4, 1);

% 舵面控制效能
P.Kr_fin          = 120;
P.Kp_fin          = 160;
P.Ky_fin          = 130;
P.B4 = [ P.Kr_fin, -P.Kr_fin,  P.Kr_fin, -P.Kr_fin;
          P.Kp_fin,  P.Kp_fin, -P.Kp_fin, -P.Kp_fin;
         -P.Ky_fin,  P.Ky_fin,  P.Ky_fin, -P.Ky_fin];
P.Ba              = P.B4 * P.Ta;
P.Ba_pinv         = pinv(P.Ba);

% 外环控制增益
P.k_alpha_q       = 7.5;
P.k_phi_p         = 4.8;
P.k_phi_d         = 2.0;
P.k_beta_r        = 7.5;
P.k_V_ud          = deg2rad(0.22);        % 速度误差到增阻指令的增益，单位 rad/(m/s)
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
P.tau_ud          = 0.20;

% 限幅参数
P.ua_max          = deg2rad([18; 20; 20]);
P.ud_max          = deg2rad(12);
P.delta_max       = deg2rad(25);
P.alpha_lim       = deg2rad([-4; 14]);
P.beta_lim        = deg2rad([-6; 6]);
P.theta_lim       = deg2rad([-25; 25]);
P.phi_lim         = deg2rad([-85; 85]);
P.gamma_lim       = deg2rad([-30; 8]);
P.V_min           = 80;
P.V_max           = 320;

% 恒值测试指令
P.n_n_cmd_const   = 0.95;
P.n_l_cmd_const   = 0.0;
P.V_cmd_const     = 190;

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

qbar_ref = 0.5 * P.rho * P.V_cmd_const^2;
CL_ref = P.n_n_cmd_const * P.m * P.g / (qbar_ref * P.S);
alpha_ref = saturate((CL_ref - P.CL0) / P.CL_alpha, P.alpha_lim(1), P.alpha_lim(2));
gamma_ref = -acos(saturate(P.n_n_cmd_const, -1.0, 1.0));
CD_base_ref = P.CD0 + P.K_induced * CL_ref^2;
CD_req_ref = P.m * P.g * max(-sin(gamma_ref), 0.0) / (qbar_ref * P.S);
P.ud_trim = saturate((CD_req_ref - CD_base_ref) / max(P.CD_ud, 1e-6), 0.0, P.ud_max);

X.z(1)            = 8000;
X.V(1)            = 200;
X.gamma(1)        = gamma_ref;
X.chi(1)          = 0.0;
X.phi(1)          = 0.0;
X.alpha(1)        = alpha_ref;
X.theta(1)        = X.gamma(1) + X.alpha(1);
X.beta(1)         = 0.0;
X.psi(1)          = X.chi(1);

%% -------------------- 结果存储 --------------------
R = struct();
R.n_n_cmd         = zeros(1, P.N);
R.n_l_cmd         = zeros(1, P.N);
R.n_total_cmd     = zeros(1, P.N);
R.V_cmd           = zeros(1, P.N);
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
R.ud_cmd          = zeros(1, P.N);
R.ud              = zeros(1, P.N);
R.delta           = zeros(4, P.N);
R.M_ctrl          = zeros(3, P.N);
R.M_aero          = zeros(3, P.N);
R.CL              = zeros(1, P.N);
R.CD              = zeros(1, P.N);
R.CD_ctrl         = zeros(1, P.N);
R.CY              = zeros(1, P.N);
R.D               = zeros(1, P.N);

[R.n_n_cmd(1), R.n_l_cmd(1), R.V_cmd(1)] = commandProfile(P.t(1), P);
R.n_total_cmd(1) = sqrt(max(R.n_n_cmd(1), 0)^2 + R.n_l_cmd(1)^2);

qbar0 = 0.5 * P.rho * X.V(1)^2;
nn_alpha_gain0 = qbar0 * P.S * P.CL_alpha / (P.m * P.g);
R.alpha_cmd(1) = saturate(R.n_total_cmd(1) / max(nn_alpha_gain0, 1e-6), ...
    P.alpha_cmd_lim(1), P.alpha_cmd_lim(2));
R.beta_cmd(1)  = 0.0;
R.phi_cmd(1)   = saturate(atan2(R.n_l_cmd(1), max(R.n_n_cmd(1), P.small_load)), ...
    -P.phi_cmd_max, P.phi_cmd_max);
R.ud_cmd(1)    = P.ud_trim;
R.ud(1)        = P.ud_trim;

[R.n_n(1), R.n_l(1), R.n_total(1), R.CL(1), R.CD(1), R.CY(1), ...
    R.CD_ctrl(1), R.D(1)] = ...
    computeLoadsAndCoeffs(X.alpha(1), X.beta(1), X.phi(1), X.V(1), P.ud_trim, P);

%% -------------------- 主仿真循环 --------------------
ua_prev = zeros(3, 1);
ud_prev = P.ud_trim;
omega_prev = [X.p(1); X.q(1); X.r(1)];
omega_dot_f_prev = zeros(3, 1);
omega_cmd_prev = zeros(3, 1);

for k = 1:P.N-1
    t = P.t(k);

    xk = [X.x(k); X.y(k); X.z(k); X.V(k); X.gamma(k); X.chi(k); ...
          X.phi(k); X.theta(k); X.psi(k); X.p(k); X.q(k); X.r(k); ...
          X.alpha(k); X.beta(k)];

    V     = X.V(k);
    alpha = X.alpha(k);
    beta  = X.beta(k);
    omega = [X.p(k); X.q(k); X.r(k)];

    qbar = 0.5 * P.rho * V^2;
    nn_alpha_gain = qbar * P.S * P.CL_alpha / (P.m * P.g);

    [n_n_cmd, n_l_cmd, V_cmd] = commandProfile(t, P);
    n_total_cmd = sqrt(max(n_n_cmd, 0)^2 + n_l_cmd^2);

    alpha_cmd = n_total_cmd / max(nn_alpha_gain, 1e-6);
    alpha_cmd = saturate(alpha_cmd, P.alpha_cmd_lim(1), P.alpha_cmd_lim(2));

    phi_cmd = atan2(n_l_cmd, max(n_n_cmd, P.small_load));
    phi_cmd = saturate(phi_cmd, -P.phi_cmd_max, P.phi_cmd_max);

    beta_cmd = 0.0;

    p_cmd = P.k_phi_p * (phi_cmd - X.phi(k)) - P.k_phi_d * X.p(k);
    q_cmd = P.k_alpha_q * (alpha_cmd - alpha);
    r_cmd = P.k_beta_r  * (beta_cmd  - beta);

    p_cmd = saturate(p_cmd, -P.p_cmd_max, P.p_cmd_max);
    q_cmd = saturate(q_cmd, -P.q_cmd_max, P.q_cmd_max);
    r_cmd = saturate(r_cmd, -P.r_cmd_max, P.r_cmd_max);
    omega_cmd = [p_cmd; q_cmd; r_cmd];

    % 根据当前滑翔工况在线计算阻力前馈，并叠加速度误差反馈
    CL_now = P.CL0 + P.CL_alpha * alpha;
    CD_base_now = P.CD0 + P.K_induced * CL_now^2 + P.CD_beta * beta^2;
    CD_req_now = P.m * P.g * max(-sin(X.gamma(k)), 0.0) / max(qbar * P.S, 1e-6);
    ud_ff = saturate((CD_req_now - CD_base_now) / max(P.CD_ud, 1e-6), 0.0, P.ud_max);
    ud_ref = saturate(ud_ff + P.k_V_ud * (V - V_cmd), 0.0, P.ud_max);
    ud_cmd = ud_prev + (P.dt / max(P.tau_ud, P.dt)) * (ud_ref - ud_prev);

    omega_cmd_dot = (omega_cmd - omega_cmd_prev) / max(P.tau_cmd_dot, P.dt);
    nu = omega_cmd_dot + P.Komega * (omega_cmd - omega);

    omega_dot_raw = (omega - omega_prev) / P.dt;
    omega_dot_meas = omega_dot_f_prev + ...
        (P.dt / P.tau_omegadot) * (omega_dot_raw - omega_dot_f_prev);

    Delta_ua = P.Ba_pinv * (P.J * (nu - omega_dot_f_prev));
    ua_cmd = ua_prev + Delta_ua;
    ua_cmd = min(max(ua_cmd, -P.ua_max), P.ua_max);

    [delta_cmd, ua_eff, ud_eff] = mixAndLimitXtail(ua_cmd, ud_cmd, P);
    M_ctrl = P.B4 * delta_cmd;

    uk.ua = ua_eff;
    uk.ud = ud_eff;
    uk.M_ctrl = M_ctrl;

    xk1 = rk4Step(@(x) plantDynamics(x, uk, P), xk, P.dt);

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

    [n_n, n_l, n_total, CL, CD, CY, CD_ctrl, D] = ...
        computeLoadsAndCoeffs(X.alpha(k+1), X.beta(k+1), X.phi(k+1), X.V(k+1), ud_eff, P);
    M_aero = aeroMoments([X.p(k); X.q(k); X.r(k)], X.alpha(k), X.beta(k), P);

    R.n_n_cmd(k)          = n_n_cmd;
    R.n_l_cmd(k)          = n_l_cmd;
    R.n_total_cmd(k)      = n_total_cmd;
    R.V_cmd(k)            = V_cmd;
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
    R.ud_cmd(k)           = ud_cmd;
    R.ud(k)               = ud_eff;
    R.delta(:, k)         = delta_cmd;
    R.M_ctrl(:, k)        = M_ctrl;
    R.M_aero(:, k)        = M_aero;
    R.CL(k+1)             = CL;
    R.CD(k+1)             = CD;
    R.CD_ctrl(k+1)        = CD_ctrl;
    R.CY(k+1)             = CY;
    R.D(k+1)              = D;

    ua_prev = ua_eff;
    ud_prev = ud_eff;
    omega_prev = omega;
    omega_dot_f_prev = omega_dot_meas;
    omega_cmd_prev = omega_cmd;
end

% 末端补齐
[R.n_n_cmd(end), R.n_l_cmd(end), R.V_cmd(end)] = commandProfile(P.t(end), P);
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
R.ud_cmd(end) = R.ud_cmd(end-1);
R.ud(end) = R.ud(end-1);
R.delta(:, end) = R.delta(:, end-1);
R.M_ctrl(:, end) = R.M_ctrl(:, end-1);
R.M_aero(:, end) = R.M_aero(:, end-1);

%% -------------------- 绘图 --------------------
t = P.t;
deg = 180 / pi;

figure(1); clf;
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(X.x / 1000, X.z / 1000, 'b', 'LineWidth', 1.8);
grid on;
xlabel('x (km)');
ylabel('z (km)');
title('Longitudinal trajectory');
nexttile;
plot3(X.x / 1000, X.y / 1000, X.z / 1000, 'r', 'LineWidth', 1.8);
grid on;
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
title('3D trajectory');
view(35, 24);

figure(2); clf;
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t, R.V_cmd, 'k--', 'LineWidth', 1.4); hold on;
plot(t, X.V, 'b', 'LineWidth', 1.8);
grid on;
xlabel('Time (s)');
ylabel('V (m/s)');
legend('V_{cmd}', 'V', 'Location', 'best');
title('Speed tracking');
nexttile;
plot(t, X.gamma * deg, 'm', 'LineWidth', 1.8);
grid on;
xlabel('Time (s)');
ylabel('\gamma (deg)');
title('Flight-path angle');

figure(3); clf;
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t, X.phi * deg, 'b', 'LineWidth', 1.8); hold on;
plot(t, R.phi_cmd * deg, 'k--', 'LineWidth', 1.3);
grid on;
xlabel('Time (s)');
ylabel('\phi (deg)');
legend('\phi', '\phi_{cmd}', 'Location', 'best');
title('Roll angle');
nexttile;
plot(t, X.theta * deg, 'b', 'LineWidth', 1.8);
grid on;
xlabel('Time (s)');
ylabel('\theta (deg)');
title('Pitch angle');
nexttile;
plot(t, X.psi * deg, 'b', 'LineWidth', 1.8);
grid on;
xlabel('Time (s)');
ylabel('\psi (deg)');
title('Yaw angle');

figure(4); clf;
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t, R.n_n_cmd, 'k--', 'LineWidth', 1.4); hold on;
plot(t, R.n_n, 'b', 'LineWidth', 1.8);
grid on;
xlabel('Time (s)');
ylabel('n_n (g)');
legend('n_{n,cmd}', 'n_n', 'Location', 'best');
title('Normal-load tracking');
nexttile;
plot(t, R.n_l_cmd, 'k--', 'LineWidth', 1.4); hold on;
plot(t, R.n_l, 'r', 'LineWidth', 1.8);
grid on;
xlabel('Time (s)');
ylabel('n_l (g)');
legend('n_{l,cmd}', 'n_l', 'Location', 'best');
title('Lateral-load tracking');
nexttile;
plot(t, R.n_total_cmd, 'k--', 'LineWidth', 1.4); hold on;
plot(t, R.n_total, 'm', 'LineWidth', 1.8);
grid on;
xlabel('Time (s)');
ylabel('n_{tot} (g)');
legend('n_{tot,cmd}', 'n_{tot}', 'Location', 'best');
title('Total load');

figure(5); clf;
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t, R.alpha_cmd * deg, 'k--', 'LineWidth', 1.4); hold on;
plot(t, X.alpha * deg, 'b', 'LineWidth', 1.8);
grid on;
xlabel('Time (s)');
ylabel('\alpha (deg)');
legend('\alpha_{cmd}', '\alpha', 'Location', 'best');
title('Angle-of-attack tracking');
nexttile;
plot(t, R.beta_cmd * deg, 'k--', 'LineWidth', 1.4); hold on;
plot(t, X.beta * deg, 'r', 'LineWidth', 1.8);
grid on;
xlabel('Time (s)');
ylabel('\beta (deg)');
legend('\beta_{cmd}', '\beta', 'Location', 'best');
title('Sideslip tracking');
nexttile;
plot(t, R.phi_cmd * deg, 'k--', 'LineWidth', 1.4); hold on;
plot(t, X.phi * deg, 'b', 'LineWidth', 1.8);
grid on;
xlabel('Time (s)');
ylabel('\phi (deg)');
legend('\phi_{cmd}', '\phi', 'Location', 'best');
title('Roll-angle tracking');

figure(6); clf;
tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:4
    nexttile;
    plot(t, R.delta(i, :) * deg, 'LineWidth', 1.8);
    grid on;
    xlabel('Time (s)');
    ylabel(sprintf('\\delta_%d (deg)', i));
    title(sprintf('Fin deflection \\delta_%d', i));
end

figure(7); clf;
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t, X.p * deg, 'b', 'LineWidth', 1.8); hold on;
plot(t, R.p_cmd * deg, 'k--', 'LineWidth', 1.3);
grid on;
xlabel('Time (s)');
ylabel('p (deg/s)');
legend('p', 'p_{cmd}', 'Location', 'best');
title('Roll rate');
nexttile;
plot(t, X.q * deg, 'b', 'LineWidth', 1.8); hold on;
plot(t, R.q_cmd * deg, 'k--', 'LineWidth', 1.3);
grid on;
xlabel('Time (s)');
ylabel('q (deg/s)');
legend('q', 'q_{cmd}', 'Location', 'best');
title('Pitch rate');
nexttile;
plot(t, X.r * deg, 'b', 'LineWidth', 1.8); hold on;
plot(t, R.r_cmd * deg, 'k--', 'LineWidth', 1.3);
grid on;
xlabel('Time (s)');
ylabel('r (deg/s)');
legend('r', 'r_{cmd}', 'Location', 'best');
title('Yaw rate');

figure(8); clf;
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t, R.ud_cmd * deg, 'k--', 'LineWidth', 1.4); hold on;
plot(t, R.ud * deg, 'b', 'LineWidth', 1.8);
grid on;
xlabel('Time (s)');
ylabel('u_d (deg)');
legend('u_{d,cmd}', 'u_d', 'Location', 'best');
title('Symmetric drag mode');
nexttile;
plot(t, R.CD, 'b', 'LineWidth', 1.8); hold on;
plot(t, R.CD_ctrl, 'r--', 'LineWidth', 1.4);
grid on;
xlabel('Time (s)');
ylabel('C_D');
legend('C_D', 'C_{D,ctrl}', 'Location', 'best');
title('Drag coefficient');
nexttile;
plot(t, R.D, 'm', 'LineWidth', 1.8);
grid on;
xlabel('Time (s)');
ylabel('D (N)');
title('Total drag');

%% -------------------- 结果摘要 --------------------
R.summary = buildSummary(R, X, P);

fprintf('\nBTT-INDI 姿态 + 阻力联合控制仿真完成。\n');
fprintf('法向过载峰值误差: %.3f g\n', R.summary.peak_nn_error_g);
fprintf('侧向过载峰值误差: %.3f g\n', R.summary.peak_nl_error_g);
fprintf('速度峰值误差: %.3f m/s\n', R.summary.peak_speed_error_mps);
fprintf('最大舵偏角: %.2f deg\n', R.summary.max_delta_deg);
fprintf('最大对称增阻模态: %.2f deg\n', R.summary.max_ud_deg);
fprintf('末端速度: %.2f m/s\n', X.V(end));
fprintf('末端俯仰角: %.2f deg\n', X.theta(end) * deg);
fprintf('末端航迹角: %.2f deg\n', X.gamma(end) * deg);
fprintf('Ba 条件数: %.2f\n', R.summary.cond_Ba);

disp(' ');
disp('Model notes:');
disp('- The BTT attitude loop is unchanged in structure: alpha controls load, phi controls load direction, beta is regulated to zero.');
disp('- A symmetric drag mode u_d is added to the X-tail mixer to support constant-speed testing.');
disp('- Attitude authority is protected first; the drag mode only uses the residual fin margin.');
disp('- alpha / beta still use low-order dynamics, so this remains a control-law verification model.');

%% ==================== 本地函数 ====================
function [n_n_cmd, n_l_cmd, V_cmd] = commandProfile(~, P)
% 恒定法向过载 + 恒定速度测试项

    n_n_cmd = P.n_n_cmd_const;
    n_l_cmd = P.n_l_cmd_const;
    V_cmd   = P.V_cmd_const;
end

function dx = plantDynamics(x, u, P)
% 简化但结构完整的平动 + 转动动力学
% x = [x; y; z; V; gamma; chi; phi; theta; psi; p; q; r; alpha; beta]

    V     = max(x(4), P.V_min);
    gamma = x(5);
    chi   = x(6);
    phi   = x(7);
    theta = x(8);
    p     = x(10);
    q     = x(11);
    r     = x(12);
    alpha = x(13);
    beta  = x(14);

    omega = [p; q; r];

    [~, ~, ~, CL, CD, CY] = computeLoadsAndCoeffs(alpha, beta, phi, V, u.ud, P); %#ok<ASGLU>
    qbar = 0.5 * P.rho * V^2;
    L = qbar * P.S * CL;
    D = qbar * P.S * CD;

    x_dot = V * cos(gamma) * cos(chi);
    y_dot = V * cos(gamma) * sin(chi);
    z_dot = V * sin(gamma);
    V_dot = -D / P.m - P.g * sin(gamma);
    gamma_dot = (L * cos(phi) / P.m - P.g * cos(gamma)) / max(V, P.V_min);
    chi_dot = (L * sin(phi) / (P.m * max(V, P.V_min) * max(cos(gamma), 0.08)));

    T = eulerRateMatrix(phi, theta);
    euler_dot = T * omega;

    M_aero = aeroMoments(omega, alpha, beta, P);
    omega_dot = P.J \ (M_aero + u.M_ctrl - cross(omega, P.J * omega));

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

function [n_n, n_l, n_total, CL, CD, CY, CD_ctrl, D] = ...
    computeLoadsAndCoeffs(alpha, beta, phi, V, ud, P)
% 含对称增阻项的 BTT 过载与气动系数模型

    qbar = 0.5 * P.rho * V^2;
    CL = P.CL0 + P.CL_alpha * alpha;
    CD_ctrl = P.CD_ud * max(ud, 0.0);
    CD = P.CD0 + P.K_induced * CL^2 + P.CD_beta * beta^2 + CD_ctrl;
    CY = P.CY_beta * beta;

    L = qbar * P.S * CL;
    D = qbar * P.S * CD;
    Y = qbar * P.S * CY;

    n_total = L / (P.m * P.g);
    n_n = n_total * cos(phi);
    n_l = n_total * sin(phi) + Y / (P.m * P.g);
end

function [delta_limited, ua_eff, ud_eff] = mixAndLimitXtail(ua_cmd, ud_cmd, P)
% 姿态控制优先，对称增阻模态仅使用剩余的正向舵偏裕度

    delta_att = P.Ta * ua_cmd;
    max_abs = max(abs(delta_att));

    if max_abs > P.delta_max
        delta_att = delta_att * (P.delta_max / max_abs);
    end

    ud_margin = min(P.delta_max - delta_att);
    ud_eff = saturate(ud_cmd, 0.0, min(P.ud_max, max(ud_margin, 0.0)));

    delta_limited = delta_att + P.Td * ud_eff;
    ua_eff = pinv(P.Ta) * delta_att;
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
% 经典四阶 Runge-Kutta 单步积分

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
% 仿真结果摘要指标

    summary = struct();
    summary.peak_nn_error_g      = max(abs(R.n_n_cmd - R.n_n));
    summary.peak_nl_error_g      = max(abs(R.n_l_cmd - R.n_l));
    summary.peak_speed_error_mps = max(abs(R.V_cmd - X.V));
    summary.rms_nn_error_g       = sqrt(mean((R.n_n_cmd - R.n_n).^2));
    summary.rms_nl_error_g       = sqrt(mean((R.n_l_cmd - R.n_l).^2));
    summary.rms_speed_error_mps  = sqrt(mean((R.V_cmd - X.V).^2));
    summary.max_delta_deg        = max(abs(R.delta(:))) * 180 / pi;
    summary.max_modal_deg        = max(abs(R.ua(:))) * 180 / pi;
    summary.max_ud_deg           = max(abs(R.ud)) * 180 / pi;
    summary.max_CD_ctrl          = max(R.CD_ctrl);
    summary.final_speed          = X.V(end);
    summary.final_range          = sqrt(X.x(end)^2 + X.y(end)^2);
    summary.cond_Ba              = cond(P.Ba);
end
