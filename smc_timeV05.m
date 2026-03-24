clc; clear; close all;

%% =========================================
% 1. 参数设置
% ==========================================
g = 9.81;
gvec = [0; 0; -g];

% ---------- 初始条件 ----------
x0 = 0;
y0 = 0;
z0 = 8000;

V0 = 250;
psi0 = deg2rad(15);
gamma0 = deg2rad(-8);

Vx0 = V0*cos(gamma0)*cos(psi0);
Vy0 = V0*cos(gamma0)*sin(psi0);
Vz0 = V0*sin(gamma0);

% ---------- 目标位置 ----------
xT = 20000;
yT = 5000;
zT = 0;
rT = [xT; yT; zT];

% ---------- 期望末端角度 ----------
psi_d = deg2rad(20);
gamma_d = deg2rad(-45);

ev_d = [cos(gamma_d)*cos(psi_d);
        cos(gamma_d)*sin(psi_d);
        sin(gamma_d)];

% ---------- 虚拟目标点参数 ----------
Ld = 3000;
rv = rT - Ld * ev_d;

% ---------- 末段切换距离 ----------
Rs = 5000;

% ---------- 基准参考速度 ----------
Vref0 = V0;

% ---------- 时间约束参数（时间误差 -> 期望速度） ----------
Td = 120;
k_v1 = 0.8;
k_v2 = 1.5;
phi_t = 10.0;

Vref_min = 180;
Vref_max = 280;

% ---------- 时间估计参数 ----------
Veff_min = 30;
beta_L = -0.05;
p_L    = 2;
k_tau  = 0.15;

% ---------- 斜率校正参数 ----------
k_slope  = 1.0;
phi_s    = 0.12;
sigma_s1 = 0.00;
sigma_s2 = 0.20;

% ---------- 第一滑模面参数（方向/角度） ----------
Lambda = diag([0.03, 0.03, 0.03]);
K1     = diag([1.6, 1.6, 1.6]);
Eta1   = diag([4.0, 4.0, 4.0]);
phi1   = [15; 15; 15];

% ---------- 第二滑模面参数（速度模长） ----------
c_v      = 0.05;     % 速度误差积分权重
k_vmag1  = 0.8;      % 速度滑模线性项增益
k_vmag2  = 1.2;      % 速度滑模切换项增益
phi_v2   = 5.0;      % 速度滑模边界层（单位：m/s）
u_t_max  = 4.0;      % 切向控制最大标量加速度

% ---------- 速度期望变化率限幅 ----------
Vdes_dot_max = 8.0;  % 期望速度变化率限幅

% ---------- 简化阻力参数 ----------
cD = 1.5e-4;

% ---------- 控制限幅 ----------
u_max = [20; 20; 20];

% ---------- 扰动 ----------
dist_amp = [0.3; 0.3; 0.2];

% ---------- 仿真参数 ----------
dt = 0.01;
tEnd = 200;
N = floor(tEnd/dt) + 1;

% ---------- 命中判据 ----------
R_hit = 20;

%% =========================================
% 2. 变量初始化
% ==========================================
X = zeros(6, N);
U = zeros(3, N);
U_raw_hist = zeros(3, N);
U_n_hist = zeros(3, N);
U_t_hist = zeros(3, N);

S1 = zeros(3, N);          % 第一滑模面
S2_hist = zeros(1, N);     % 第二滑模面（速度模长）
zv_hist = zeros(1, N);     % 速度误差积分状态

R_hist = zeros(1, N);
sigma_hist = zeros(1, N);
ang_err_hist = zeros(1, N);
psi_hist = zeros(1, N);
gamma_hist = zeros(1, N);

% ---------- 视线角记录 ----------
psi_los_hist = zeros(1, N);
gamma_los_hist = zeros(1, N);

% ---------- 时间约束记录 ----------
tgo_hist = zeros(1, N);
tgo_raw_hist = zeros(1, N);
tgo_dot_hist = zeros(1, N);
et_hist = zeros(1, N);
Vr_hist = zeros(1, N);
Veff_hist = zeros(1, N);
Lgo_hist = zeros(1, N);
kc_hist = zeros(1, N);
Vdes_hist = zeros(1, N);
Vdes_dot_hist = zeros(1, N);

% ---------- 时间通道拆分记录 ----------
Rt_hist = zeros(1, N);
Lextra_hist = zeros(1, N);
Veff_t_hist = zeros(1, N);
kc_t_hist = zeros(1, N);
dtL_hist = zeros(1, N);
dtchi_hist = zeros(1, N);
dts_hist = zeros(1, N);
t_rv_hist = zeros(1, N);
misalign_hist = zeros(1, N);
mode_hist = zeros(1, N);
ws_hist = zeros(1, N);
slope_err_hist = zeros(1, N);

% ---------- 速度模长记录 ----------
Vnorm_hist = zeros(1, N);

% ---------- 三轴速度分量记录 ----------
Vx_hist = zeros(1, N);
Vy_hist = zeros(1, N);
Vz_hist = zeros(1, N);

% ---------- 惯性系加速度分量记录 ----------
A_cmd_hist = zeros(3, N);
A_drag_hist = zeros(3, N);
A_total_hist = zeros(3, N);

% ---------- 速度系下加速度分量记录 ----------
A_cmd_v_hist = zeros(3, N);
A_total_v_hist = zeros(3, N);

% ---------- 第二滑模面相关记录 ----------
eVmag_hist = zeros(1, N);

% ---------- 动态期望速度状态 ----------
Vdes_state = Vref0;

% ---------- 误差记录 ----------
er_hist = zeros(3, N);
eV_hist = zeros(3, N);

X(:,1) = [x0; y0; z0; Vx0; Vy0; Vz0];

% ---------- 速度误差积分状态 ----------
z_v = 0;

hit_flag = false;
hit_index = N;
ground_hit_flag = false;

%% =========================================
% 3. 主循环
% ==========================================
for k = 1:N-1

    % 当前时间
    t = (k-1) * dt;

    % 当前状态
    r = X(1:3, k);
    V = X(4:6, k);
    Vnorm = norm(V);
    if Vnorm < 1e-6
        Vnorm = 1e-6;
    end
    Vnorm_hist(k) = Vnorm;

    % ---------- 记录三轴速度分量 ----------
    Vx_hist(k) = V(1);
    Vy_hist(k) = V(2);
    Vz_hist(k) = V(3);

    % 当前速度方向
    ev_cur = V / Vnorm;

    % 相对位置
    Rvec = rT - r;
    R = norm(Rvec);
    R_hist(k) = R;
    Rt_hist(k) = R;

    % 当前视线方向
    if R < 1e-6
        eR = [1;0;0];
    else
        eR = Rvec / R;
    end

    % ---------- 记录视线角 ----------
    psi_los_hist(k) = atan2(Rvec(2), Rvec(1));
    gamma_los_hist(k) = atan2(Rvec(3), sqrt(Rvec(1)^2 + Rvec(2)^2));

    % ---------- 末段切换因子 ----------
    if R >= Rs
        sigma = 0;
    else
        sigma = 1 - R / Rs;
    end
    sigma = max(0, min(1, sigma));
    sigma_hist(k) = sigma;

    % =========================================
    % A. 角度/位置参考方向：混合参考点
    % =========================================
    r_ref = (1 - sigma) * rv + sigma * rT;
    ref_vec = r_ref - r;
    ref_norm = max(norm(ref_vec), 1e-6);
    eref = ref_vec / ref_norm;

    % =========================================
    % B. 时间通道：R/Vr + Δt_L + Δt_χ + Δt_s
    % =========================================
    Vr = dot(V, eR);
    Vr_hist(k) = Vr;

    Veff_t = max(Vr, Veff_min);
    Veff_t_hist(k) = Veff_t;
    Veff_hist(k) = Veff_t;

    % ---------- 1) 主体项 ----------
    t_rv = R / Veff_t;
    t_rv_hist(k) = t_rv;

    % ---------- 2) 附加路径修正 ----------
    Lextra = beta_L * (1 - sigma)^p_L * Ld;
    Lextra_hist(k) = Lextra;

    V_L = max(Vnorm, Veff_min);
    dt_L = Lextra / V_L;
    dtL_hist(k) = dt_L;

    % ---------- 3) 曲率修正 ----------
    misalign = 1 - dot(ev_cur, eR);
    misalign = max(0, misalign);
    misalign_hist(k) = misalign;

    dt_chi = k_tau * misalign * R / max(Vnorm, Veff_min);
    dtchi_hist(k) = dt_chi;

    % ---------- 原始剩余时间 ----------
    tgo_raw = t_rv + dt_L + dt_chi;
    tgo_raw_hist(k) = tgo_raw;

    if k == 1
        tgo_dot = -1;
    else
        tgo_dot = (tgo_raw - tgo_raw_hist(k-1)) / dt;
    end
    tgo_dot_hist(k) = tgo_dot;

    % ---------- 切换附近权重 ----------
    if sigma <= sigma_s1
        ws = 0;
    elseif sigma >= sigma_s2
        ws = 1;
    else
        ws = (sigma - sigma_s1) / (sigma_s2 - sigma_s1);
    end
    ws_hist(k) = ws;

    slope_err = tgo_dot + 1;
    slope_err_hist(k) = slope_err;

    dt_s = -ws * k_slope * slope_err;
    dt_s = max(min(dt_s, phi_s), -phi_s);
    dts_hist(k) = dt_s;

    % ---------- 最终剩余时间 ----------
    tgo = tgo_raw + dt_s;
    tgo_hist(k) = tgo;

    % ---------- 兼容记录 ----------
    Lgo = R + Lextra;
    Lgo_hist(k) = Lgo;

    kc_eq = tgo / max(t_rv, 1e-6);
    kc_t_hist(k) = kc_eq;
    kc_hist(k) = kc_eq;
    mode_hist(k) = 0;

    % ---------- 到达时间误差 ----------
    et = t + tgo - Td;
    et_hist(k) = et;

    % =========================================
    % C. 时间误差 -> 期望速度
    % =========================================
    sat_et = sat_func(et / phi_t);
    Vdes_dot = k_v1 * et + k_v2 * sat_et;
    Vdes_dot = max(min(Vdes_dot, Vdes_dot_max), -Vdes_dot_max);
    Vdes_dot_hist(k) = Vdes_dot;

    Vdes_state = Vdes_state + Vdes_dot * dt;
    Vdes_state = min(max(Vdes_state, Vref_min), Vref_max);
    Vdes = Vdes_state;
    Vdes_hist(k) = Vdes;

    % =========================================
    % D. 第一滑模面：方向/角度
    % =========================================
    % 只用期望方向构造参考速度，速度模长取当前速度模长
    Vd_dir = Vnorm * eref;

    if k == 1
        Vd_dir_dot = [0; 0; 0];
    else
        r_prev = X(1:3, k-1);
        R_prev = norm(rT - r_prev);

        if R_prev >= Rs
            sigma_prev = 0;
        else
            sigma_prev = 1 - R_prev / Rs;
        end
        sigma_prev = max(0, min(1, sigma_prev));

        r_ref_prev = (1 - sigma_prev) * rv + sigma_prev * rT;
        ref_vec_prev = r_ref_prev - r_prev;
        eref_prev = ref_vec_prev / max(norm(ref_vec_prev), 1e-6);

        Vnorm_prev = max(norm(X(4:6,k-1)), 1e-6);
        Vd_dir_prev = Vnorm_prev * eref_prev;
        Vd_dir_dot = (Vd_dir - Vd_dir_prev) / dt;
    end

    er = r - rT;
    eV = V - Vd_dir;
    er_hist(:,k) = er;
    eV_hist(:,k) = eV;

    s1 = eV + Lambda * er;
    S1(:,k) = s1;

    sat_s1 = zeros(3,1);
    for i = 1:3
        sat_s1(i) = sat_func(s1(i)/phi1(i));
    end

    u_raw = Vd_dir_dot ...
            - Lambda * V ...
            - gvec ...
            - K1 * s1 ...
            - Eta1 * sat_s1;
    U_raw_hist(:,k) = u_raw;

    % ---------- 法向投影控制 ----------
    P_n = eye(3) - ev_cur * ev_cur.';
    u_n = P_n * u_raw;
    U_n_hist(:,k) = u_n;

    % =========================================
    % E. 第二滑模面：速度模长
    % =========================================
    eVmag = Vnorm - Vdes;
    eVmag_hist(k) = eVmag;

    z_v = z_v + eVmag * dt;
    zv_hist(k) = z_v;

    s2 = eVmag + c_v * z_v;
    S2_hist(k) = s2;

    sat_s2 = sat_func(s2 / phi_v2);

    % 当速度高于期望时，u_t 应沿 -ev_cur 减速；低于期望时沿 +ev_cur 加速
    a_t_scalar = -(k_vmag1 * s2 + k_vmag2 * sat_s2);
    a_t_scalar = max(min(a_t_scalar, u_t_max), -u_t_max);

    u_t = a_t_scalar * ev_cur;
    U_t_hist(:,k) = u_t;

    % =========================================
    % F. 合成控制
    % =========================================
    u = u_n + u_t;

    % 控制限幅
    u = max(min(u, u_max), -u_max);
    U(:,k) = u;

    % ---------- 扰动 ----------
    d = [dist_amp(1)*sin(0.20*t);
         dist_amp(2)*cos(0.15*t);
         dist_amp(3)*sin(0.35*t)];

    % ---------- 简化阻力 ----------
    a_drag = -cD * Vnorm * V;

    % ---------- 总加速度（惯性系） ----------
    a_total = u + gvec + d + a_drag;

    % ---------- 记录惯性系加速度分量 ----------
    A_cmd_hist(:,k) = u;
    A_drag_hist(:,k) = a_drag;
    A_total_hist(:,k) = a_total;

    % =========================================
    % 速度系构造
    % =========================================
    ex_v = ev_cur;

    ezI = [0; 0; 1];
    ey_v = cross(ezI, ex_v);

    if norm(ey_v) < 1e-6
        ey_v = [0; 1; 0];
    else
        ey_v = ey_v / norm(ey_v);
    end

    ez_v = cross(ex_v, ey_v);
    ez_v = ez_v / max(norm(ez_v), 1e-6);

    C_iv = [ex_v, ey_v, ez_v];

    % ---------- 速度系下加速度分量 ----------
    a_cmd_v = C_iv' * u;
    a_total_v = C_iv' * a_total;

    A_cmd_v_hist(:,k) = a_cmd_v;
    A_total_v_hist(:,k) = a_total_v;

    % ---------- 欧拉积分 ----------
    r_dot = V;
    V_dot = a_total;
    X_next = X(:,k) + [r_dot; V_dot] * dt;

    % ---------- 记录速度角 ----------
    psi_hist(k) = atan2(V(2), V(1));
    gamma_hist(k) = atan2(V(3), sqrt(V(1)^2 + V(2)^2));

    cos_ang = dot(ev_cur, ev_d);
    cos_ang = min(max(cos_ang, -1), 1);
    ang_err_hist(k) = rad2deg(acos(cos_ang));

    % ---------- 命中判断 ----------
    if R <= R_hit
        hit_flag = true;
        X(:,k+1) = X_next;
        hit_index = k+1;
        break;
    end

    % ---------- 地面交点判定 ----------
    if X(3,k) > 0 && X_next(3) <= 0
        alpha = X(3,k) / (X(3,k) - X_next(3));
        X_ground = X(:,k) + alpha * (X_next - X(:,k));
        X_ground(3) = 0;
        X(:,k+1) = X_ground;
        hit_index = k+1;
        ground_hit_flag = true;
        break;
    end

    X(:,k+1) = X_next;
end

%% =========================================
% 4. 末端数据计算与历史量补齐
% ==========================================
r_final = X(1:3, hit_index);
V_final = X(4:6, hit_index);

Vf = norm(V_final);
if Vf < 1e-6
    Vf = 1e-6;
end
ev_final = V_final / Vf;
Vnorm_hist(hit_index) = Vf;

% ---------- 补齐末端速度分量 ----------
Vx_hist(hit_index) = V_final(1);
Vy_hist(hit_index) = V_final(2);
Vz_hist(hit_index) = V_final(3);

R_final = norm(rT - r_final);
R_hist(hit_index) = R_final;
Rt_hist(hit_index) = R_final;

Rvec_final = rT - r_final;
psi_los_hist(hit_index) = atan2(Rvec_final(2), Rvec_final(1));
gamma_los_hist(hit_index) = atan2(Rvec_final(3), sqrt(Rvec_final(1)^2 + Rvec_final(2)^2));

if R_final >= Rs
    sigma_final = 0;
else
    sigma_final = 1 - R_final / Rs;
end
sigma_final = max(0, min(1, sigma_final));
sigma_hist(hit_index) = sigma_final;

% ---------- 角度/位置通道参考点 ----------
r_ref_final = (1 - sigma_final) * rv + sigma_final * rT;
ref_vec_final = r_ref_final - r_final;
ref_norm_final = max(norm(ref_vec_final), 1e-6);
eref_final = ref_vec_final / ref_norm_final;

% ---------- 真实目标方向 ----------
if R_final < 1e-6
    eR_final = [1;0;0];
else
    eR_final = Rvec_final / R_final;
end

Vr_final = dot(V_final, eR_final);
Vr_hist(hit_index) = Vr_final;

Veff_t_final = max(Vr_final, Veff_min);
Veff_t_hist(hit_index) = Veff_t_final;
Veff_hist(hit_index) = Veff_t_final;

% ---------- 1) 主体项：真实目标直达时间 ----------
t_rv_final = R_final / Veff_t_final;
t_rv_hist(hit_index) = t_rv_final;

% ---------- 2) 附加路径修正 ----------
Lextra_final = beta_L * (1 - sigma_final)^p_L * Ld;
Lextra_hist(hit_index) = Lextra_final;

V_L_final = max(Vf, Veff_min);
dt_L_final = Lextra_final / V_L_final;
dtL_hist(hit_index) = dt_L_final;

% ---------- 3) 曲率修正 ----------
misalign_final = 1 - dot(ev_final, eR_final);
misalign_final = max(0, misalign_final);
misalign_hist(hit_index) = misalign_final;

dt_chi_final = k_tau * misalign_final * R_final / max(Vf, Veff_min);
dtchi_hist(hit_index) = dt_chi_final;

% ---------- 原始剩余时间 ----------
tgo_raw_final = t_rv_final + dt_L_final + dt_chi_final;
tgo_raw_hist(hit_index) = tgo_raw_final;

if hit_index == 1
    tgo_dot_final = -1;
else
    tgo_dot_final = (tgo_raw_final - tgo_raw_hist(max(hit_index-1,1))) / dt;
end
tgo_dot_hist(hit_index) = tgo_dot_final;

if sigma_final <= sigma_s1
    ws_final = 0;
elseif sigma_final >= sigma_s2
    ws_final = 1;
else
    ws_final = (sigma_final - sigma_s1) / (sigma_s2 - sigma_s1);
end
ws_hist(hit_index) = ws_final;

slope_err_final = tgo_dot_final + 1;
slope_err_hist(hit_index) = slope_err_final;

dt_s_final = -ws_final * k_slope * slope_err_final;
dt_s_final = max(min(dt_s_final, phi_s), -phi_s);
dts_hist(hit_index) = dt_s_final;

% ---------- 最终剩余时间 ----------
tgo_final = tgo_raw_final + dt_s_final;
tgo_hist(hit_index) = tgo_final;

% ---------- 兼容记录 ----------
Lgo_final = R_final + Lextra_final;
Lgo_hist(hit_index) = Lgo_final;

kc_t_final = tgo_final / max(t_rv_final, 1e-6);
kc_t_hist(hit_index) = kc_t_final;
kc_hist(hit_index) = kc_t_final;

mode_hist(hit_index) = 0;

t_final = (hit_index - 1) * dt;
et_final = t_final + tgo_final - Td;
et_hist(hit_index) = et_final;

Vdes_hist(hit_index) = Vdes_hist(max(hit_index-1,1));
Vdes_dot_hist(hit_index) = Vdes_dot_hist(max(hit_index-1,1));

psi_f = atan2(V_final(2), V_final(1));
gamma_f = atan2(V_final(3), sqrt(V_final(1)^2 + V_final(2)^2));
psi_hist(hit_index) = psi_f;
gamma_hist(hit_index) = gamma_f;

cos_ang_f = dot(ev_final, ev_d);
cos_ang_f = min(max(cos_ang_f, -1), 1);
ang_err_f = rad2deg(acos(cos_ang_f));
ang_err_hist(hit_index) = ang_err_f;

Vd_final = Vdes_hist(hit_index) * eref_final;
er_final = r_final - rT;
eV_final = V_final - Vd_final;
s_final = eV_final + Lambda * er_final;

er_hist(:, hit_index) = er_final;
eV_hist(:, hit_index) = eV_final;
S1(:, hit_index) = s_final;

eVmag_hist(hit_index) = Vf - Vdes_hist(hit_index);
zv_hist(hit_index) = zv_hist(max(hit_index-1,1));
S2_hist(hit_index) = eVmag_hist(hit_index) + c_v * zv_hist(hit_index);

if hit_index >= 2
    U(:, hit_index) = U(:, hit_index-1);
    U_raw_hist(:, hit_index) = U_raw_hist(:, hit_index-1);
    U_n_hist(:, hit_index) = U_n_hist(:, hit_index-1);
    U_t_hist(:, hit_index) = U_t_hist(:, hit_index-1);

    A_cmd_hist(:, hit_index) = A_cmd_hist(:, hit_index-1);
    A_drag_hist(:, hit_index) = A_drag_hist(:, hit_index-1);
    A_total_hist(:, hit_index) = A_total_hist(:, hit_index-1);
    A_cmd_v_hist(:, hit_index) = A_cmd_v_hist(:, hit_index-1);
    A_total_v_hist(:, hit_index) = A_total_v_hist(:, hit_index-1);
end

%% =========================================
% 5. 截取有效数据
% ==========================================
idx = 1:hit_index;
time = (idx-1) * dt;

X = X(:, idx);
U = U(:, idx);
U_raw_hist = U_raw_hist(:, idx);
U_n_hist = U_n_hist(:, idx);
U_t_hist = U_t_hist(:, idx);
S1 = S1(:, idx);
S2_hist = S2_hist(idx);

R_hist = R_hist(idx);
sigma_hist = sigma_hist(idx);
ang_err_hist = ang_err_hist(idx);
psi_hist = psi_hist(idx);
gamma_hist = gamma_hist(idx);

psi_los_hist = psi_los_hist(idx);
gamma_los_hist = gamma_los_hist(idx);

tgo_hist = tgo_hist(idx);
tgo_raw_hist = tgo_raw_hist(idx);
tgo_dot_hist = tgo_dot_hist(idx);
et_hist = et_hist(idx);
Vr_hist = Vr_hist(idx);
Veff_hist = Veff_hist(idx);
Lgo_hist = Lgo_hist(idx);
kc_hist = kc_hist(idx);
Vdes_hist = Vdes_hist(idx);
Vdes_dot_hist = Vdes_dot_hist(idx);

Rt_hist = Rt_hist(idx);
Lextra_hist = Lextra_hist(idx);
Veff_t_hist = Veff_t_hist(idx);
kc_t_hist = kc_t_hist(idx);
dtL_hist = dtL_hist(idx);
dtchi_hist = dtchi_hist(idx);
dts_hist = dts_hist(idx);
t_rv_hist = t_rv_hist(idx);
misalign_hist = misalign_hist(idx);
mode_hist = mode_hist(idx);
ws_hist = ws_hist(idx);
slope_err_hist = slope_err_hist(idx);

Vnorm_hist = Vnorm_hist(idx);

Vx_hist = Vx_hist(idx);
Vy_hist = Vy_hist(idx);
Vz_hist = Vz_hist(idx);

A_cmd_hist = A_cmd_hist(:, idx);
A_drag_hist = A_drag_hist(:, idx);
A_total_hist = A_total_hist(:, idx);
A_cmd_v_hist = A_cmd_v_hist(:, idx);
A_total_v_hist = A_total_v_hist(:, idx);

er_hist = er_hist(:, idx);
eV_hist = eV_hist(:, idx);
eVmag_hist = eVmag_hist(idx);
zv_hist = zv_hist(idx);

r_final = X(1:3,end);
V_final = X(4:6,end);
Vf = norm(V_final);
if Vf < 1e-6
    Vf = 1e-6;
end
ev_final = V_final / Vf;

psi_f = atan2(V_final(2), V_final(1));
gamma_f = atan2(V_final(3), sqrt(V_final(1)^2 + V_final(2)^2));

%% =========================================
% 6. 脱靶量
% ==========================================
dr = rT - r_final;
lambda_star = dot(dr, ev_final);

if lambda_star >= 0
    miss_distance = norm(dr - lambda_star * ev_final);
else
    miss_distance = norm(dr);
end

cos_ang_f = dot(ev_final, ev_d);
cos_ang_f = min(max(cos_ang_f, -1), 1);
ang_err_f = rad2deg(acos(cos_ang_f));

actual_tf = time(end);
time_error = actual_tf - Td;

%% =========================================
% 7. 输出结果
% ==========================================
fprintf('================ 仿真结果 ================\n');
fprintf('是否命中: %d\n', hit_flag);
fprintf('是否地面截获: %d\n', ground_hit_flag);
fprintf('末端位置: [%.2f, %.2f, %.2f] m\n', r_final(1), r_final(2), r_final(3));
fprintf('目标位置: [%.2f, %.2f, %.2f] m\n', xT, yT, zT);
fprintf('脱靶量: %.4f m\n', miss_distance);
fprintf('末端速度: %.4f m/s\n', norm(V_final));
fprintf('末端航向角 psi_f   = %.4f deg\n', rad2deg(psi_f));
fprintf('期望航向角 psi_d   = %.4f deg\n', rad2deg(psi_d));
fprintf('末端弹道倾角 gamma_f = %.4f deg\n', rad2deg(gamma_f));
fprintf('期望弹道倾角 gamma_d = %.4f deg\n', rad2deg(gamma_d));
fprintf('末端速度方向误差角 = %.4f deg\n', ang_err_f);
fprintf('指定到达时间 Td   = %.4f s\n', Td);
fprintf('实际到达时间 tf   = %.4f s\n', actual_tf);
fprintf('到达时间误差 tf-Td = %.4f s\n', time_error);
fprintf('末端期望速度 Vdes = %.4f m/s\n', Vdes_hist(end));

%% =========================================
% 8. 绘图
% ==========================================

figure;
plot3(X(1,:), X(2,:), X(3,:), 'b', 'LineWidth', 1.8); hold on;
plot3(xT, yT, zT, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
plot3(rv(1), rv(2), rv(3), 'ks', 'MarkerSize', 7, 'LineWidth', 1.5);
grid on;
xlabel('x / m'); ylabel('y / m'); zlabel('z / m');
title('三维飞行轨迹');
legend('飞行轨迹', '目标', '虚拟目标');
set(gca, 'ZDir', 'normal');

figure;
plot(time, Vnorm_hist, 'LineWidth', 1.5); hold on;
plot(time, Vdes_hist, '--r', 'LineWidth', 1.2);
grid on;
xlabel('t / s'); ylabel('V / m/s');
title('速度模长与期望速度');
legend('|V|', 'Vdes');

figure;
subplot(6,1,1);
plot(time, tgo_hist, 'LineWidth', 1.5); hold on;
plot(time, tgo_raw_hist, '--', 'LineWidth', 1.1);
grid on;
ylabel('t_{go}');
title('剩余时间估计：R/V_r + \Deltat_L + \Deltat_\chi + \Deltat_s');
legend('tgo', 'tgo\_raw');

subplot(6,1,2);
plot(time, et_hist, 'LineWidth', 1.5); hold on;
yline(0, '--r', 'LineWidth', 1.0);
grid on;
ylabel('e_t');

subplot(6,1,3);
plot(time, t_rv_hist, 'LineWidth', 1.5); hold on;
plot(time, dtL_hist, '--', 'LineWidth', 1.2);
plot(time, dtchi_hist, ':', 'LineWidth', 1.2);
plot(time, dts_hist, '-.', 'LineWidth', 1.2);
grid on;
ylabel('时间 / s');
legend('R/V_r', '\Deltat_L', '\Deltat_\chi', '\Deltat_s');

subplot(6,1,4);
plot(time, tgo_dot_hist, 'LineWidth', 1.5); hold on;
yline(-1, '--r', 'LineWidth', 1.0);
grid on;
ylabel('dot(t_{go})');

subplot(6,1,5);
plot(time, ws_hist, 'LineWidth', 1.5); hold on;
plot(time, slope_err_hist, '--', 'LineWidth', 1.2);
grid on;
ylabel('w_s / err');
legend('w_s', 'slope\_err');

subplot(6,1,6);
plot(time, misalign_hist, 'LineWidth', 1.5); grid on;
ylabel('1-e_v^Te_R');
xlabel('t / s');

figure;
subplot(3,1,1);
plot(time, Vnorm_hist - Vdes_hist, 'LineWidth', 1.5); hold on;
yline(0,'--r');
grid on;
ylabel('e_{Vmag}');
title('速度模长误差');

subplot(3,1,2);
plot(time, zv_hist, 'LineWidth', 1.5); grid on;
ylabel('z_v');
title('速度误差积分状态');

subplot(3,1,3);
plot(time, S2_hist, 'LineWidth', 1.5); hold on;
yline(0,'--r');
grid on;
ylabel('s_2');
xlabel('t / s');
title('第二滑模面（速度模长）');

figure;
subplot(2,1,1);
plot(time, rad2deg(psi_los_hist), 'LineWidth', 1.5); grid on;
ylabel('\psi_{LOS} / deg');
title('视线方位角变化');

subplot(2,1,2);
plot(time, rad2deg(gamma_los_hist), 'LineWidth', 1.5); grid on;
ylabel('\gamma_{LOS} / deg');
xlabel('t / s');
title('视线俯仰角变化');

figure;
subplot(3,1,1);
plot(time, Vx_hist, 'LineWidth', 1.5); grid on;
ylabel('V_x / m/s');
title('三轴速度分量');

subplot(3,1,2);
plot(time, Vy_hist, 'LineWidth', 1.5); grid on;
ylabel('V_y / m/s');

subplot(3,1,3);
plot(time, Vz_hist, 'LineWidth', 1.5); grid on;
ylabel('V_z / m/s');
xlabel('t / s');

figure;
subplot(3,1,1);
plot(time, A_cmd_hist(1,:), 'LineWidth', 1.5); grid on;
ylabel('a_{cx} / m/s^2');
title('惯性系下控制加速度分量');

subplot(3,1,2);
plot(time, A_cmd_hist(2,:), 'LineWidth', 1.5); grid on;
ylabel('a_{cy} / m/s^2');

subplot(3,1,3);
plot(time, A_cmd_hist(3,:), 'LineWidth', 1.5); grid on;
ylabel('a_{cz} / m/s^2');
xlabel('t / s');

figure;
subplot(3,1,1);
plot(time, A_drag_hist(1,:), 'LineWidth', 1.5); grid on;
ylabel('a_{dx} / m/s^2');
title('惯性系下阻力加速度分量');

subplot(3,1,2);
plot(time, A_drag_hist(2,:), 'LineWidth', 1.5); grid on;
ylabel('a_{dy} / m/s^2');

subplot(3,1,3);
plot(time, A_drag_hist(3,:), 'LineWidth', 1.5); grid on;
ylabel('a_{dz} / m/s^2');
xlabel('t / s');

figure;
subplot(3,1,1);
plot(time, A_total_hist(1,:), 'LineWidth', 1.5); grid on;
ylabel('a_x / m/s^2');
title('惯性系下总加速度分量');

subplot(3,1,2);
plot(time, A_total_hist(2,:), 'LineWidth', 1.5); grid on;
ylabel('a_y / m/s^2');

subplot(3,1,3);
plot(time, A_total_hist(3,:), 'LineWidth', 1.5); grid on;
ylabel('a_z / m/s^2');
xlabel('t / s');

figure;
subplot(3,1,1);
plot(time, A_cmd_v_hist(1,:), 'LineWidth', 1.5); grid on;
ylabel('a_{cvx} / m/s^2');
title('速度系下控制加速度分量');

subplot(3,1,2);
plot(time, A_cmd_v_hist(2,:), 'LineWidth', 1.5); grid on;
ylabel('a_{cvy} / m/s^2');

subplot(3,1,3);
plot(time, A_cmd_v_hist(3,:), 'LineWidth', 1.5); grid on;
ylabel('a_{cvz} / m/s^2');
xlabel('t / s');

figure;
subplot(3,1,1);
plot(time, A_total_v_hist(1,:), 'LineWidth', 1.5); grid on;
ylabel('a_{vx} / m/s^2');
title('速度系下总加速度分量');

subplot(3,1,2);
plot(time, A_total_v_hist(2,:), 'LineWidth', 1.5); grid on;
ylabel('a_{vy} / m/s^2');

subplot(3,1,3);
plot(time, A_total_v_hist(3,:), 'LineWidth', 1.5); grid on;
ylabel('a_{vz} / m/s^2');
xlabel('t / s');

%% =========================================
% 9. 饱和函数
% ==========================================
function y = sat_func(x)
    if x > 1
        y = 1;
    elseif x < -1
        y = -1;
    else
        y = x;
    end
end