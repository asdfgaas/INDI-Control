%% glider_time_indi_sim.m
% Unpowered glider control and simulation under strict time-of-arrival constraint
%
% Model structure:
% 1) 3D point-mass translational model in velocity coordinates.
% 2) Simplified attitude loop with first-order theta/phi tracking.
% 3) BTT logic is reserved in the lateral channel through phi_c, but the
%    main verification in this script focuses on the longitudinal channel.
%
% Control logic:
% 1) Time-of-arrival is the highest-priority objective.
% 2) Normal acceleration tracking is the second-priority objective.
% 3) Speed reference is only an auxiliary shaping target and is suppressed
%    when the time error is large.
% 4) Longitudinal control uses a single dominant input alpha_c and solves
%    an INDI-style weighted least-squares increment problem based on online
%    estimates of g_n = da_n/du and g_t = da_t/du.
%
% Important note:
% This script intentionally stops at the attitude/rate-loop vicinity.
% It computes theta_c / phi_c and equivalent q_c / p_c, but does not model
% rigid-body rotational dynamics, actuator dynamics, or X-tail mixing.

clear; clc; close all;

%% -------------------- Parameters --------------------
P = struct();

% Environment and vehicle
P.g            = 9.81;          % m/s^2
P.m            = 900;           % kg
P.S            = 1.60;          % m^2
P.rho0         = 1.225;         % kg/m^3
P.h_scale      = 8500;          % m

% Aerodynamics
P.CL0          = 0.10;
P.CLa          = 4.40;          % 1/rad
P.CD0          = 0.038;
P.K            = 0.065;

% Inner-loop approximations
P.tau_theta    = 0.28;          % s
P.tau_phi      = 0.25;          % s

% Limits and numerical protections
P.alpha_min    = deg2rad(-3);
P.alpha_max    = deg2rad(11);
P.theta_min    = deg2rad(-30);
P.theta_max    = deg2rad( 30);
P.phi_min      = deg2rad(-50);
P.phi_max      = deg2rad( 50);
P.phi_cmd_lim  = deg2rad(20);
P.du_max       = deg2rad(0.60); % max alpha increment per control step
P.V_min        = 90;            % m/s
P.cosGamMin    = 0.05;
P.cosPhiMin    = 0.20;
P.range_hit    = 80;            % m
P.z_ground     = 0;

% Guidance-like references
P.a_n_c        = 2.5;           % m/s^2, constant normal acceleration command
P.a_t_c        = 0.0;           % m/s^2, feedforward tangential acceleration command
P.V_ref        = 210;           % m/s, auxiliary speed reference only
P.T_arrival    = 100;           % s, desired arrival time

% Time-dominant tangential channel
P.k_tau        = 0.16;          % main gain for time error
P.k_tau_i      = 0.035;         % integral gain for residual time error
P.k_vg_tau     = 0.06;          % desired ground-speed shaping from time-to-go
P.k_V_t        = 0.045;         % auxiliary speed shaping gain
P.tau_gate     = 5.0;           % s, gate width for sigma_t
P.e_tau_int_lim = 60.0;         % s, anti-windup bound for integrated time error

% INDI tracking gains and weights
P.k_n_indi     = 0.70;
P.k_t_indi     = 0.95;
P.w_t_hi       = 5.0;           % time channel weight when |e_tau| is large
P.w_t_lo       = 1.6;           % time channel weight when |e_tau| is small
P.w_n          = 1.0;           % normal acceleration weight
P.lambda_u     = 0.22;          % control regularization

% Online estimation / filtering
P.tau_y_f      = 0.18;          % s, first-order filter for a_n and a_t
P.rls_lambda   = 0.985;         % forgetting factor
P.P0           = 25.0;          % initial covariance for scalar RLS
P.du_update_th = deg2rad(0.08); % do not update if input motion is too small
P.gn_min       = 15.0;          % lower bound of da_n/du
P.gn_max       = 180.0;         % upper bound of da_n/du
P.gt_min       = -20.0;         % lower bound of da_t/du
P.gt_max       = -0.20;         % upper bound of da_t/du

% BTT lateral loop (kept simple here)
P.psi_ref      = 0.0;           % rad
P.k_psi        = 1.30;          % psi error to phi_c

% Target point and simulation setup
P.r_target     = [16500; 0; 8000]; % m, timing target projected in horizontal mission plane
P.dt_ctrl      = 0.05;          % s
P.t_end        = 110;           % s

% Tuned notes:
% This parameter set was adjusted after end-to-end MATLAB runs to keep the
% main scenario stable and to better show the difference between time-driven
% energy management and the comparison case.

%% -------------------- Initial Conditions --------------------
z0             = 8000;          % m
V0             = 235;           % m/s
gamma0         = deg2rad(-5.0); % rad
psi0           = 0.0;           % rad
phi0           = 0.0;           % rad
alpha_trim     = deg2rad(2.5);  % rad
theta0         = gamma0 + alpha_trim;
x0             = 0.0;
y0             = 0.0;

X0 = [x0; y0; z0; V0; gamma0; psi0; theta0; phi0];

%% -------------------- Scenario Definitions --------------------
SC(1).name = 'Time-dominant energy management';
SC(1).enable_time_priority = true;
SC(1).color = [0.00, 0.45, 0.74];
SC(1).line_style = '-';

SC(2).name = 'No time-dominant term';
SC(2).enable_time_priority = false;
SC(2).color = [0.85, 0.33, 0.10];
SC(2).line_style = '--';

%% -------------------- Simulation --------------------
RES = cell(numel(SC), 1);
for i = 1:numel(SC)
    RES{i} = runScenario(P, X0, SC(i));
end

%% -------------------- Plotting --------------------
plotResults(RES, P, SC);

%% -------------------- Summary --------------------
fprintf('\n========== Time-Constrained Glider Simulation Summary ==========\n');
for i = 1:numel(SC)
    printSummary(RES{i}, P, SC(i));
end

%% ==================== Local Functions ====================
function R = runScenario(P, X0, SC)
% Discrete-time controller + segmented ode45 propagation

    t_vec = 0:P.dt_ctrl:P.t_end;
    N = numel(t_vec);

    X = zeros(8, N);
    X(:, 1) = X0;

    % Histories
    alpha_c_hist = zeros(1, N);
    theta_c_hist = zeros(1, N);
    phi_c_hist   = zeros(1, N);
    q_c_hist     = zeros(1, N);
    p_c_hist     = zeros(1, N);
    r_c_hist     = zeros(1, N);

    a_n_hist     = zeros(1, N);
    a_t_hist     = zeros(1, N);
    a_n_f_hist   = zeros(1, N);
    a_t_f_hist   = zeros(1, N);
    a_t_star_hist = zeros(1, N);
    t_go_hist    = zeros(1, N);
    t_d_hist     = zeros(1, N);
    e_tau_hist   = zeros(1, N);
    sigma_t_hist = zeros(1, N);
    e_n_hist     = zeros(1, N);
    e_t_hist     = zeros(1, N);
    e_tau_int_hist = zeros(1, N);
    g_n_hist     = zeros(1, N);
    g_t_hist     = zeros(1, N);
    range_hist   = zeros(1, N);
    Vg_des_hist  = zeros(1, N);

    % Initial outputs
    S0 = outputModel(X0, P);
    a_n_hist(1)  = S0.a_n;
    a_t_hist(1)  = S0.a_t;
    a_n_f_hist(1) = S0.a_n;
    a_t_f_hist(1) = S0.a_t;

    % Initial estimator guesses from local aerodynamic slopes
    est.g_n = saturate(S0.gn_model, P.gn_min, P.gn_max);
    est.g_t = saturate(S0.gt_model, P.gt_min, P.gt_max);
    est.Pn  = P.P0;
    est.Pt  = P.P0;

    alpha_cmd_prev = saturate(X0(7) - X0(5), P.alpha_min, P.alpha_max);
    alpha_c_hist(1) = alpha_cmd_prev;
    theta_c_hist(1) = X0(7);
    phi_c_hist(1)   = X0(8);
    g_n_hist(1)     = est.g_n;
    g_t_hist(1)     = est.g_t;

    range_hist(1)   = norm(P.r_target(1:2) - X0(1:2));
    t_go_hist(1)    = range_hist(1) / max(X0(4) * cos(X0(5)), 0.55 * P.V_min);
    t_d_hist(1)     = P.T_arrival;
    e_tau_hist(1)   = t_go_hist(1) - t_d_hist(1);
    sigma_t_hist(1) = saturate(1.0 - abs(e_tau_hist(1)) / P.tau_gate, 0.0, 1.0)^2;
    Vg_des_hist(1)  = range_hist(1) / max(t_d_hist(1), P.dt_ctrl);
    e_tau_int = 0.0;

    stop_idx = N;

    for k = 1:N-1
        t = t_vec(k);
        xk = X(:, k);
        Sk = outputModel(xk, P);

        % Remaining-time quantities
        pos = xk(1:3);
        % Time-to-go is defined from the horizontal remaining range here.
        % This keeps the timing problem well-posed for the current point-mass
        % model, which does not include terminal dive guidance or full 3D
        % impact-angle shaping.
        range_to_go = norm(P.r_target(1:2) - pos(1:2));
        V_eff_go = max(xk(4) * cos(xk(5)), 0.55 * P.V_min);
        t_go = range_to_go / V_eff_go;
        t_d = max(P.T_arrival - t, 0);
        e_tau = t_go - t_d;
        e_tau_int = saturate(e_tau_int + e_tau * P.dt_ctrl, -P.e_tau_int_lim, P.e_tau_int_lim);
        Vg_des = range_to_go / max(t_d, P.dt_ctrl);
        e_vg_tau = Vg_des - V_eff_go;

        % Time-priority scheduling:
        % sigma_t -> 0 when time error is large, so speed reference is suppressed.
        sigma_t = saturate(1.0 - abs(e_tau) / P.tau_gate, 0.0, 1.0);
        sigma_t = sigma_t^2;

        if SC.enable_time_priority
            a_t_star = P.a_t_c + P.k_tau * e_tau + P.k_tau_i * e_tau_int ...
                + P.k_vg_tau * e_vg_tau ...
                - P.k_V_t * sigma_t * (xk(4) - P.V_ref);
        else
            a_t_star = P.a_t_c - P.k_V_t * (xk(4) - P.V_ref);
        end

        % BTT lateral channel is kept, but the current case uses psi_ref = 0.
        psi_err = wrapToPiLocal(P.psi_ref - xk(6));
        phi_c = saturate(P.k_psi * psi_err, -P.phi_cmd_lim, P.phi_cmd_lim);

        % INDI weighted least-squares increment with one dominant input u = alpha_c
        e_n = P.a_n_c - a_n_f_hist(k);
        e_t = a_t_star - a_t_f_hist(k);
        dy_d = [P.k_t_indi * e_t; P.k_n_indi * e_n];

        w_t = P.w_t_hi * (1 - sigma_t) + P.w_t_lo * sigma_t;
        W2 = diag([w_t^2, P.w_n^2]);
        g_vec = [est.g_t; est.g_n];

        delta_u = (g_vec' * W2 * dy_d) / (g_vec' * W2 * g_vec + P.lambda_u);
        delta_u = saturate(delta_u, -P.du_max, P.du_max);

        alpha_c = saturate(alpha_cmd_prev + delta_u, P.alpha_min, P.alpha_max);
        theta_c = saturate(alpha_c + xk(5), P.theta_min, P.theta_max);

        q_c = (theta_c - xk(7)) / P.tau_theta;
        p_c = (phi_c - xk(8)) / P.tau_phi;
        r_c = 0.0;

        % Propagate the point-mass dynamics over one control interval using ode45
        U.theta_c = theta_c;
        U.phi_c   = phi_c;
        U.alpha_c = alpha_c;

        [~, x_seg] = ode45(@(~, xx) pointMassDyn(xx, U, P), [t, t_vec(k+1)], xk);
        X(:, k+1) = x_seg(end, :)';

        % New outputs and filtered outputs
        Sk1 = outputModel(X(:, k+1), P);
        a_n_hist(k+1) = Sk1.a_n;
        a_t_hist(k+1) = Sk1.a_t;
        a_n_f_hist(k+1) = a_n_f_hist(k) + (P.dt_ctrl / P.tau_y_f) * (a_n_hist(k+1) - a_n_f_hist(k));
        a_t_f_hist(k+1) = a_t_f_hist(k) + (P.dt_ctrl / P.tau_y_f) * (a_t_hist(k+1) - a_t_f_hist(k));

        % Online scalar RLS updates on output increments
        du = alpha_c - alpha_cmd_prev;
        if abs(du) > P.du_update_th
            delta_an = a_n_f_hist(k+1) - a_n_f_hist(k);
            delta_at = a_t_f_hist(k+1) - a_t_f_hist(k);

            [est.g_n, est.Pn] = rlsScalar(est.g_n, est.Pn, du, delta_an, P.rls_lambda);
            [est.g_t, est.Pt] = rlsScalar(est.g_t, est.Pt, du, delta_at, P.rls_lambda);

            est.g_n = saturate(est.g_n, P.gn_min, P.gn_max);
            est.g_t = saturate(est.g_t, P.gt_min, P.gt_max);
        end

        % Record histories
        alpha_c_hist(k+1) = alpha_c;
        theta_c_hist(k+1) = theta_c;
        phi_c_hist(k+1)   = phi_c;
        q_c_hist(k+1)     = q_c;
        p_c_hist(k+1)     = p_c;
        r_c_hist(k+1)     = r_c;
        a_t_star_hist(k+1) = a_t_star;
        t_go_hist(k+1)    = t_go;
        t_d_hist(k+1)     = t_d;
        e_tau_hist(k+1)   = e_tau;
        sigma_t_hist(k+1) = sigma_t;
        e_n_hist(k+1)     = e_n;
        e_t_hist(k+1)     = e_t;
        e_tau_int_hist(k+1) = e_tau_int;
        g_n_hist(k+1)     = est.g_n;
        g_t_hist(k+1)     = est.g_t;
        range_hist(k+1)   = range_to_go;
        Vg_des_hist(k+1)  = Vg_des;

        alpha_cmd_prev = alpha_c;

        % Stop if target is reached or the vehicle hits the ground
        range_next = norm(P.r_target(1:2) - X(1:2, k+1));
        if range_next <= P.range_hit || X(3, k+1) <= P.z_ground
            stop_idx = k + 1;
            break;
        end
    end

    idx = 1:stop_idx;
    X = X(:, idx);
    t_vec = t_vec(idx);

    R.t         = t_vec(:);
    R.x         = X(1, idx).';
    R.y         = X(2, idx).';
    R.z         = X(3, idx).';
    R.V         = X(4, idx).';
    R.gamma     = X(5, idx).';
    R.psi       = X(6, idx).';
    R.theta     = X(7, idx).';
    R.phi       = X(8, idx).';
    R.alpha     = R.theta - R.gamma;
    R.alpha_c   = alpha_c_hist(idx).';
    R.theta_c   = theta_c_hist(idx).';
    R.phi_c     = phi_c_hist(idx).';
    R.q_c       = q_c_hist(idx).';
    R.p_c       = p_c_hist(idx).';
    R.r_c       = r_c_hist(idx).';
    R.a_n       = a_n_hist(idx).';
    R.a_t       = a_t_hist(idx).';
    R.a_n_f     = a_n_f_hist(idx).';
    R.a_t_f     = a_t_f_hist(idx).';
    R.a_t_star  = a_t_star_hist(idx).';
    R.t_go      = t_go_hist(idx).';
    R.t_d       = t_d_hist(idx).';
    R.e_tau     = e_tau_hist(idx).';
    R.sigma_t   = sigma_t_hist(idx).';
    R.e_n       = e_n_hist(idx).';
    R.e_t       = e_t_hist(idx).';
    R.e_tau_int = e_tau_int_hist(idx).';
    R.g_n       = g_n_hist(idx).';
    R.g_t       = g_t_hist(idx).';
    R.range     = range_hist(idx).';
    R.Vg_des    = Vg_des_hist(idx).';
    R.arrival_time = R.t(end);
    R.arrival_error = R.arrival_time - P.T_arrival;
    R.name      = SC.name;
end

function dx = pointMassDyn(x, U, P)
% Point-mass translational dynamics + simplified fast attitude loop

    z     = x(3);
    V     = max(x(4), P.V_min);
    gamma = x(5);
    psi   = x(6);
    theta = x(7);
    phi   = x(8);

    rho = airDensity(z, P);
    alpha = theta - gamma;
    CL = P.CL0 + P.CLa * alpha;
    CD = P.CD0 + P.K * CL^2;
    qbar = 0.5 * rho * V^2;
    L = qbar * P.S * CL;
    D = qbar * P.S * CD;

    cosGamma = signedLowerBound(cos(gamma), P.cosGamMin);

    x_dot     = V * cos(gamma) * cos(psi);
    y_dot     = V * cos(gamma) * sin(psi);
    z_dot     = V * sin(gamma);
    V_dot     = -D / P.m - P.g * sin(gamma);
    gamma_dot = (L * cos(phi) / P.m - P.g * cos(gamma)) / V;
    psi_dot   = (L * sin(phi) / (P.m * V * cosGamma));

    theta_dot = (U.theta_c - theta) / P.tau_theta;
    phi_dot   = (U.phi_c   - phi  ) / P.tau_phi;

    dx = [x_dot; y_dot; z_dot; V_dot; gamma_dot; psi_dot; theta_dot; phi_dot];
end

function S = outputModel(x, P)
% Compute outputs and local incremental control effectiveness

    z     = x(3);
    V     = max(x(4), P.V_min);
    gamma = x(5);
    theta = x(7);
    phi   = x(8);

    rho = airDensity(z, P);
    alpha = theta - gamma;
    CL = P.CL0 + P.CLa * alpha;
    CD = P.CD0 + P.K * CL^2;
    qbar = 0.5 * rho * V^2;
    coeff = qbar * P.S / P.m;

    L = qbar * P.S * CL;
    D = qbar * P.S * CD;

    S.a_n = L * cos(phi) / P.m - P.g * cos(gamma);
    S.a_t = -D / P.m - P.g * sin(gamma);
    S.CL  = CL;
    S.CD  = CD;
    S.rho = rho;

    % Local slopes with u = alpha_c ~= alpha for the inner-loop-dominant model
    S.gn_model = coeff * P.CLa * cos(phi);
    S.gt_model = -coeff * (2 * P.K * CL * P.CLa);
end

function [g_new, P_new] = rlsScalar(g_old, P_old, du, dy, lambda)
% Scalar recursive least squares for dy = g * du

    denom = lambda + du * P_old * du;
    K = (P_old * du) / denom;
    g_new = g_old + K * (dy - g_old * du);
    P_new = (P_old - K * du * P_old) / lambda;
end

function rho = airDensity(z, P)
    h = max(z, 0.0);
    rho = P.rho0 * exp(-h / P.h_scale);
end

function y = saturate(u, umin, umax)
    y = min(max(u, umin), umax);
end

function y = signedLowerBound(u, boundAbs)
    if abs(u) < boundAbs
        y = sign(u + eps) * boundAbs;
    else
        y = u;
    end
end

function ang = wrapToPiLocal(ang)
    ang = mod(ang + pi, 2*pi) - pi;
end

function plotResults(RES, P, SC)
    deg = 180 / pi;
    td_plot = RES{1}.t_d;

    % Figure 1: time-of-arrival and energy-management quantities
    figure('Color', 'w', 'Name', 'Time and Energy Management');
    tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plot(RES{1}.t, RES{1}.t_go, 'Color', SC(1).color, 'LineWidth', 1.8); hold on;
    plot(RES{1}.t, td_plot, 'k:', 'LineWidth', 1.6);
    plot(RES{2}.t, RES{2}.t_go, 'Color', SC(2).color, 'LineWidth', 1.8, 'LineStyle', SC(2).line_style);
    grid on;
    xlabel('t (s)');
    ylabel('Time (s)');
    title('t_{go} and t_d');
    legend('t_{go} main', 't_d', 't_{go} compare', 'Location', 'best');

    nexttile;
    plot(RES{1}.t, RES{1}.e_tau, 'Color', SC(1).color, 'LineWidth', 1.8); hold on;
    plot(RES{2}.t, RES{2}.e_tau, 'Color', SC(2).color, 'LineWidth', 1.8, 'LineStyle', SC(2).line_style);
    yline(0, 'k:', 'LineWidth', 1.5);
    grid on;
    xlabel('t (s)');
    ylabel('e_{\tau} (s)');
    title('Time Error');
    legend('Main case', 'Compare case', 'Location', 'best');

    nexttile;
    plot(RES{1}.t, RES{1}.sigma_t, 'Color', SC(1).color, 'LineWidth', 1.8); hold on;
    plot(RES{2}.t, RES{2}.sigma_t, 'Color', SC(2).color, 'LineWidth', 1.8, 'LineStyle', SC(2).line_style);
    grid on;
    xlabel('t (s)');
    ylabel('\sigma_t');
    title('Time-Gating Factor');
    legend('Main case', 'Compare case', 'Location', 'best');

    nexttile;
    plot(RES{1}.t, RES{1}.a_n, 'Color', SC(1).color, 'LineWidth', 1.8); hold on;
    plot(RES{2}.t, RES{2}.a_n, 'Color', SC(2).color, 'LineWidth', 1.8, 'LineStyle', SC(2).line_style);
    yline(P.a_n_c, 'k:', 'LineWidth', 1.5);
    grid on;
    xlabel('t (s)');
    ylabel('a_n (m/s^2)');
    title('Normal Acceleration Tracking');
    legend('Main case', 'Compare case', 'a_{n,c}', 'Location', 'best');

    nexttile;
    plot(RES{1}.t, RES{1}.a_t, 'Color', SC(1).color, 'LineWidth', 1.8); hold on;
    plot(RES{1}.t, RES{1}.a_t_star, 'k:', 'LineWidth', 1.6);
    plot(RES{2}.t, RES{2}.a_t, 'Color', SC(2).color, 'LineWidth', 1.8, 'LineStyle', SC(2).line_style);
    grid on;
    xlabel('t (s)');
    ylabel('a_t (m/s^2)');
    title('Tangential Channel');
    legend('a_t main', 'a_t^* main', 'a_t compare', 'Location', 'best');

    sgtitle('Time Constraint and Energy Management');

    % Figure 2: longitudinal state and control response
    figure('Color', 'w', 'Name', 'Longitudinal Control Response');
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plot(RES{1}.t, RES{1}.V, 'Color', SC(1).color, 'LineWidth', 1.8); hold on;
    plot(RES{2}.t, RES{2}.V, 'Color', SC(2).color, 'LineWidth', 1.8, 'LineStyle', SC(2).line_style);
    yline(P.V_ref, 'k:', 'LineWidth', 1.5);
    grid on;
    xlabel('t (s)');
    ylabel('V (m/s)');
    title('Speed');
    legend('Main case', 'Compare case', 'V_{ref}', 'Location', 'best');

    nexttile;
    plot(RES{1}.t, RES{1}.alpha * deg, 'Color', SC(1).color, 'LineWidth', 1.8); hold on;
    plot(RES{1}.t, RES{1}.alpha_c * deg, 'k:', 'LineWidth', 1.6);
    plot(RES{2}.t, RES{2}.alpha * deg, 'Color', SC(2).color, 'LineWidth', 1.8, 'LineStyle', SC(2).line_style);
    grid on;
    xlabel('t (s)');
    ylabel('\alpha (deg)');
    title('Angle of Attack');
    legend('\alpha main', '\alpha_c main', '\alpha compare', 'Location', 'best');

    nexttile;
    plot(RES{1}.t, RES{1}.theta * deg, 'Color', SC(1).color, 'LineWidth', 1.8); hold on;
    plot(RES{1}.t, RES{1}.gamma * deg, 'Color', [0.47, 0.67, 0.19], 'LineWidth', 1.8);
    plot(RES{2}.t, RES{2}.theta * deg, 'Color', SC(2).color, 'LineWidth', 1.8, 'LineStyle', SC(2).line_style);
    grid on;
    xlabel('t (s)');
    ylabel('Angle (deg)');
    title('\theta and \gamma');
    legend('\theta main', '\gamma main', '\theta compare', 'Location', 'best');

    nexttile;
    plot(RES{1}.t, RES{1}.q_c * deg, 'Color', SC(1).color, 'LineWidth', 1.8); hold on;
    plot(RES{1}.t, RES{1}.p_c * deg, 'k--', 'LineWidth', 1.6);
    grid on;
    xlabel('t (s)');
    ylabel('Rate cmd (deg/s)');
    title('Equivalent Rate Commands');
    legend('q_c', 'p_c', 'Location', 'best');

    sgtitle('Longitudinal Control and Inner-Loop Commands');

    % Figure 3: lateral channel and trajectory
    figure('Color', 'w', 'Name', 'Attitude and Trajectory');
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plot(RES{1}.t, RES{1}.phi * deg, 'Color', SC(1).color, 'LineWidth', 1.8); hold on;
    plot(RES{1}.t, RES{1}.phi_c * deg, 'k:', 'LineWidth', 1.6);
    plot(RES{2}.t, RES{2}.phi * deg, 'Color', SC(2).color, 'LineWidth', 1.8, 'LineStyle', SC(2).line_style);
    grid on;
    xlabel('t (s)');
    ylabel('\phi (deg)');
    title('BTT Roll Channel');
    legend('\phi main', '\phi_c main', '\phi compare', 'Location', 'best');

    nexttile;
    plot(RES{1}.t, RES{1}.z, 'Color', SC(1).color, 'LineWidth', 1.8); hold on;
    plot(RES{2}.t, RES{2}.z, 'Color', SC(2).color, 'LineWidth', 1.8, 'LineStyle', SC(2).line_style);
    grid on;
    xlabel('t (s)');
    ylabel('z (m)');
    title('Altitude');
    legend('Main case', 'Compare case', 'Location', 'best');

    nexttile;
    plot(RES{1}.x / 1000, RES{1}.z / 1000, 'Color', SC(1).color, 'LineWidth', 1.8); hold on;
    plot(RES{2}.x / 1000, RES{2}.z / 1000, 'Color', SC(2).color, 'LineWidth', 1.8, 'LineStyle', SC(2).line_style);
    plot(P.r_target(1) / 1000, P.r_target(3) / 1000, 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
    grid on;
    xlabel('x (km)');
    ylabel('z (km)');
    title('x-z Trajectory');
    legend('Main case', 'Compare case', 'Target', 'Location', 'best');

    nexttile;
    plot3(RES{1}.x / 1000, RES{1}.y / 1000, RES{1}.z / 1000, 'Color', SC(1).color, 'LineWidth', 1.8); hold on;
    plot3(RES{2}.x / 1000, RES{2}.y / 1000, RES{2}.z / 1000, 'Color', SC(2).color, 'LineWidth', 1.8, 'LineStyle', SC(2).line_style);
    plot3(P.r_target(1) / 1000, P.r_target(2) / 1000, P.r_target(3) / 1000, 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
    grid on;
    xlabel('x (km)');
    ylabel('y (km)');
    zlabel('z (km)');
    title('3D Trajectory');
    legend('Main case', 'Compare case', 'Target', 'Location', 'best');
    view(35, 25);

    sgtitle('Attitude Channel and Flight Path');

    % Figure 4: estimator diagnostics
    figure('Color', 'w', 'Name', 'INDI Estimation Diagnostics');
    tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plot(RES{1}.t, RES{1}.g_n, 'Color', SC(1).color, 'LineWidth', 1.8); hold on;
    plot(RES{1}.t, RES{1}.g_t, 'k--', 'LineWidth', 1.6);
    grid on;
    xlabel('t (s)');
    ylabel('Estimate');
    title('Online INDI Effectiveness Estimates');
    legend('g_n', 'g_t', 'Location', 'best');

    nexttile;
    plot(RES{1}.t, RES{1}.e_tau_int, 'Color', SC(1).color, 'LineWidth', 1.8); hold on;
    plot(RES{1}.t, RES{1}.e_t, 'k--', 'LineWidth', 1.6);
    grid on;
    xlabel('t (s)');
    ylabel('Diagnostic');
    title('Time-Channel Internal Signals');
    legend('Integrated e_{\tau}', 'e_t = a_t^* - a_t', 'Location', 'best');

    sgtitle('INDI Estimation and Time-Channel Diagnostics');
end

function printSummary(R, P, SC)
    fprintf('Scenario: %s\n', SC.name);
    fprintf('  Final time             : %7.2f s\n', R.t(end));
    fprintf('  Estimated t_go final   : %7.2f s\n', R.t_go(end));
    fprintf('  Arrival time error     : %+7.2f s (actual - desired)\n', R.arrival_error);
    fprintf('  Final time error e_tau : %+7.2f s\n', R.e_tau(end));
    fprintf('  Final a_n              : %7.2f m/s^2 (cmd %5.2f)\n', R.a_n(end), P.a_n_c);
    fprintf('  Final V                : %7.2f m/s   (ref %6.2f)\n', R.V(end), P.V_ref);
    fprintf('  Final altitude         : %7.2f m\n', R.z(end));
    fprintf('  Final horizontal range : %7.2f m\n', norm(P.r_target(1:2) - [R.x(end); R.y(end)]));
    fprintf('  Max |alpha|            : %7.2f deg\n', max(abs(R.alpha)) * 180 / pi);
    fprintf('  Max |q_c|              : %7.2f deg/s\n\n', max(abs(R.q_c)) * 180 / pi);
end
