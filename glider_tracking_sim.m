%% glider_tracking_sim.m
% Unpowered glider point-mass tracking simulation
%
% Modeling notes:
% 1) A 3D point-mass glider model is used. Only translational motion and
%    the states V, gamma, psi, theta, phi, and position are retained.
% 2) No rigid-body rotational dynamics are modeled. The script does not
%    integrate p/q/r dynamics and does not compute control surface deflection.
% 3) The control structure follows a simplified BTT-style idea:
%    - normal acceleration command -> baseline alpha command
%    - speed command -> alpha correction term
%    - alpha command -> theta command through alpha = theta - gamma
% 4) The inner loop is approximated by first-order attitude tracking:
%       theta_dot = (theta_c - theta)/tau_theta
%       phi_dot   = (phi_c   - phi  )/tau_phi
%    This is equivalent to designing only to the attitude/rate-loop vicinity.
%
% Stability note:
% The parameter set in this script was tuned to remain robust for the main
% case. If a_n_ref is increased to 8~10 m/s^2 under unpowered flight,
% speed loss becomes much stronger and alpha saturation may occur earlier.

clear; clc; close all;

%% -------------------- Parameter Definition --------------------
P = struct();

% Constants and vehicle data
P.g          = 9.81;            % m/s^2
P.m          = 900;             % kg
P.S          = 1.60;            % m^2
P.rho0       = 1.225;           % kg/m^3
P.h_scale    = 8500;            % m, exponential atmosphere scale height
P.V_min      = 80;              % m/s, numerical floor for V
P.cosGamMin  = 0.05;            % avoid division by zero in psi_dot

% Aerodynamic data (alpha in rad)
P.CL0        = 0.10;
P.CLa        = 4.50;            % 1/rad
P.CD0        = 0.04;
P.K          = 0.070;

% Control parameters
P.tau_theta  = 0.35;            % s, ideal pitch inner-loop time constant
P.tau_phi    = 0.30;            % s, ideal roll inner-loop time constant
P.kV         = 0.0025;          % rad/(m/s), speed error to alpha correction
P.phi_c      = 0.0;             % rad, reserved roll command

% Limits
P.alpha_min  = deg2rad(-3);     % rad
P.alpha_max  = deg2rad(10);     % rad
P.theta_min  = deg2rad(-30);    % rad
P.theta_max  = deg2rad( 30);    % rad
P.phi_min    = deg2rad(-60);    % rad
P.phi_max    = deg2rad( 60);    % rad

% Commands
P.V_ref      = 195;             % m/s
P.a_n_ref    = 3.5;             % m/s^2

% Simulation horizon
t_span       = [0, 100];         % s

%% -------------------- Initial Condition --------------------
z0           = 8000;            % m, altitude (positive upward)
V0           = 220;             % m/s
gamma0       = deg2rad(-5);     % rad
psi0         = 0.0;             % rad
phi0         = 0.0;             % rad

alpha_trim   = deg2rad(3.5);    % rad
theta0       = gamma0 + alpha_trim;

x0           = 0.0;
y0           = 0.0;

% State vector: [x; y; z; V; gamma; psi; theta; phi]
X0 = [x0; y0; z0; V0; gamma0; psi0; theta0; phi0];

%% -------------------- Scenario Definition --------------------
% Scenario 1: normal acceleration + speed correction
SC1 = struct();
SC1.name            = 'With speed correction';
SC1.enable_speed_fb = true;
SC1.color           = [0.00, 0.45, 0.74];
SC1.line_style      = '-';

% Scenario 2: only normal acceleration tracking, no speed correction
SC2 = struct();
SC2.name            = 'Without speed correction';
SC2.enable_speed_fb = false;
SC2.color           = [0.85, 0.33, 0.10];
SC2.line_style      = '--';

%% -------------------- Numerical Integration --------------------
ode_opts = odeset( ...
    'RelTol', 1e-7, ...
    'AbsTol', 1e-8, ...
    'Events', @(t, X) stopEvents(t, X, P));

[t1, X1] = ode45(@(t, X) pointMassModel(t, X, P, SC1), t_span, X0, ode_opts);
[t2, X2] = ode45(@(t, X) pointMassModel(t, X, P, SC2), t_span, X0, ode_opts);

R1 = postProcess(t1, X1, P, SC1);
R2 = postProcess(t2, X2, P, SC2);

%% -------------------- Plotting --------------------
deg = 180 / pi;

figure('Color', 'w', 'Name', 'Glider Tracking Simulation');
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 1. Speed tracking
nexttile;
plot(R1.t, R1.V, 'Color', SC1.color, 'LineWidth', 1.8); hold on;
plot(R2.t, R2.V, 'Color', SC2.color, 'LineWidth', 1.8, 'LineStyle', SC2.line_style);
plot(R1.t, R1.V_ref * ones(size(R1.t)), 'k:', 'LineWidth', 1.6);
grid on;
xlabel('t (s)');
ylabel('V (m/s)');
title('Speed Tracking');
legend('Main case', 'No speed feedback', 'V_{ref}', 'Location', 'best');

% 2. Normal acceleration tracking
nexttile;
plot(R1.t, R1.a_n, 'Color', SC1.color, 'LineWidth', 1.8); hold on;
plot(R2.t, R2.a_n, 'Color', SC2.color, 'LineWidth', 1.8, 'LineStyle', SC2.line_style);
plot(R1.t, R1.a_n_ref * ones(size(R1.t)), 'k:', 'LineWidth', 1.6);
grid on;
xlabel('t (s)');
ylabel('a_n (m/s^2)');
title('Normal Acceleration Tracking');
legend('Main case', 'No speed feedback', 'a_{n,ref}', 'Location', 'best');

% 3. Alpha and alpha command
nexttile;
plot(R1.t, R1.alpha * deg, 'Color', SC1.color, 'LineWidth', 1.8); hold on;
plot(R1.t, R1.alpha_c * deg, 'Color', [0.10, 0.10, 0.10], 'LineWidth', 1.5, 'LineStyle', ':');
plot(R2.t, R2.alpha * deg, 'Color', SC2.color, 'LineWidth', 1.8, 'LineStyle', SC2.line_style);
grid on;
xlabel('t (s)');
ylabel('\alpha (deg)');
title('Angle of Attack');
legend('\alpha (main)', '\alpha_c (main)', '\alpha (no speed fb)', 'Location', 'best');

% 4. Theta and gamma
nexttile;
plot(R1.t, R1.theta * deg, 'Color', SC1.color, 'LineWidth', 1.8); hold on;
plot(R1.t, R1.gamma * deg, 'Color', [0.47, 0.67, 0.19], 'LineWidth', 1.8);
plot(R2.t, R2.theta * deg, 'Color', SC2.color, 'LineWidth', 1.6, 'LineStyle', SC2.line_style);
grid on;
xlabel('t (s)');
ylabel('Angle (deg)');
title('\theta and \gamma');
legend('\theta (main)', '\gamma (main)', '\theta (no speed fb)', 'Location', 'best');

% 5. Altitude
nexttile;
plot(R1.t, R1.z, 'Color', SC1.color, 'LineWidth', 1.8); hold on;
plot(R2.t, R2.z, 'Color', SC2.color, 'LineWidth', 1.8, 'LineStyle', SC2.line_style);
grid on;
xlabel('t (s)');
ylabel('z (m)');
title('Altitude');
legend('Main case', 'No speed feedback', 'Location', 'best');

% 6. x-z trajectory
nexttile;
plot(R1.x / 1000, R1.z / 1000, 'Color', SC1.color, 'LineWidth', 1.8); hold on;
plot(R2.x / 1000, R2.z / 1000, 'Color', SC2.color, 'LineWidth', 1.8, 'LineStyle', SC2.line_style);
grid on;
xlabel('x (km)');
ylabel('z (km)');
title('Longitudinal Trajectory');
legend('Main case', 'No speed feedback', 'Location', 'best');

sgtitle('Unpowered Glider Tracking: Constant V_{ref} and a_{n,ref}');

% Optional extra plots for phi and 3D path
figure('Color', 'w', 'Name', 'Additional Plots');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(R1.t, R1.phi * deg, 'Color', SC1.color, 'LineWidth', 1.8); hold on;
plot(R1.t, R1.phi_c * deg, 'k:', 'LineWidth', 1.5);
plot(R2.t, R2.phi * deg, 'Color', SC2.color, 'LineWidth', 1.8, 'LineStyle', SC2.line_style);
grid on;
xlabel('t (s)');
ylabel('\phi (deg)');
title('Roll Channel');
legend('\phi (main)', '\phi_c', '\phi (no speed fb)', 'Location', 'best');

nexttile;
plot3(R1.x / 1000, R1.y / 1000, R1.z / 1000, 'Color', SC1.color, 'LineWidth', 1.8); hold on;
plot3(R2.x / 1000, R2.y / 1000, R2.z / 1000, 'Color', SC2.color, 'LineWidth', 1.8, 'LineStyle', SC2.line_style);
grid on;
axis equal;
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
title('3D Trajectory');
legend('Main case', 'No speed feedback', 'Location', 'best');
view(35, 25);

%% -------------------- Summary Print --------------------
fprintf('\n========== Simulation Summary ==========\n');
printSummary(R1, P, SC1);
printSummary(R2, P, SC2);

%% ==================== Local Functions ====================
function dX = pointMassModel(~, X, P, SC)
% 3D point-mass model with simplified attitude inner loop

    z     = X(3);
    V     = X(4);
    gamma = X(5);
    psi   = X(6);
    theta = X(7);
    phi   = X(8);

    rho = airDensity(z, P);
    Veff = max(V, P.V_min);
    cosGamma = signedLowerBound(cos(gamma), P.cosGamMin);

    % Current aerodynamics
    alpha = theta - gamma;
    CL = P.CL0 + P.CLa * alpha;
    CD = P.CD0 + P.K * CL^2;
    qbar = 0.5 * rho * Veff^2;
    L = qbar * P.S * CL;
    D = qbar * P.S * CD;

    % Actual normal acceleration as required in the problem statement
    a_n = L * cos(phi) / P.m - P.g * cos(gamma);

    % Baseline alpha from normal acceleration command
    cosPhi = signedLowerBound(cos(phi), 0.2);
    L_cmd = P.m * (P.a_n_ref + P.g * cos(gamma)) / cosPhi;
    CL_cmd = 2 * L_cmd / max(rho * Veff^2 * P.S, 1.0);
    alpha_n = (CL_cmd - P.CL0) / P.CLa;

    % Speed correction: if V is higher than V_ref, increase alpha to add drag;
    % if V is lower than V_ref, decrease alpha to relieve drag.
    if SC.enable_speed_fb
        delta_alpha_V = P.kV * (V - P.V_ref);
    else
        delta_alpha_V = 0.0;
    end

    alpha_c = saturate(alpha_n + delta_alpha_V, P.alpha_min, P.alpha_max);
    theta_c = saturate(alpha_c + gamma, P.theta_min, P.theta_max);
    phi_c = saturate(P.phi_c, P.phi_min, P.phi_max);

    % Ideal attitude inner loop
    theta_dot = (theta_c - theta) / P.tau_theta;
    phi_dot   = (phi_c   - phi  ) / P.tau_phi;

    % Translational equations
    x_dot     = Veff * cos(gamma) * cos(psi);
    y_dot     = Veff * cos(gamma) * sin(psi);
    z_dot     = Veff * sin(gamma);
    V_dot     = -D / P.m - P.g * sin(gamma);
    gamma_dot = a_n / Veff;
    psi_dot   = (L * sin(phi) / (P.m * Veff * cosGamma));

    dX = [x_dot; y_dot; z_dot; V_dot; gamma_dot; psi_dot; theta_dot; phi_dot];
end

function R = postProcess(t, X, P, SC)
% Compute control and aerodynamic variables for plots and analysis

    n = numel(t);

    R.t        = t;
    R.x        = X(:, 1);
    R.y        = X(:, 2);
    R.z        = X(:, 3);
    R.V        = X(:, 4);
    R.gamma    = X(:, 5);
    R.psi      = X(:, 6);
    R.theta    = X(:, 7);
    R.phi      = X(:, 8);
    R.alpha    = zeros(n, 1);
    R.alpha_c  = zeros(n, 1);
    R.alpha_n  = zeros(n, 1);
    R.theta_c  = zeros(n, 1);
    R.phi_c    = zeros(n, 1);
    R.a_n      = zeros(n, 1);
    R.CL       = zeros(n, 1);
    R.CD       = zeros(n, 1);
    R.rho      = zeros(n, 1);
    R.q_c      = zeros(n, 1);
    R.p_c      = zeros(n, 1);
    R.r_c      = zeros(n, 1);
    R.V_ref    = P.V_ref;
    R.a_n_ref  = P.a_n_ref;
    R.name     = SC.name;

    for k = 1:n
        z     = X(k, 3);
        V     = max(X(k, 4), P.V_min);
        gamma = X(k, 5);
        theta = X(k, 7);
        phi   = X(k, 8);

        rho = airDensity(z, P);
        alpha = theta - gamma;
        CL = P.CL0 + P.CLa * alpha;
        CD = P.CD0 + P.K * CL^2;
        qbar = 0.5 * rho * V^2;
        L = qbar * P.S * CL;

        a_n = L * cos(phi) / P.m - P.g * cos(gamma);

        cosPhi = signedLowerBound(cos(phi), 0.2);
        L_cmd = P.m * (P.a_n_ref + P.g * cos(gamma)) / cosPhi;
        CL_cmd = 2 * L_cmd / max(rho * V^2 * P.S, 1.0);
        alpha_n = (CL_cmd - P.CL0) / P.CLa;

        if SC.enable_speed_fb
            delta_alpha_V = P.kV * (V - P.V_ref);
        else
            delta_alpha_V = 0.0;
        end

        alpha_c = saturate(alpha_n + delta_alpha_V, P.alpha_min, P.alpha_max);
        theta_c = saturate(alpha_c + gamma, P.theta_min, P.theta_max);
        phi_c = saturate(P.phi_c, P.phi_min, P.phi_max);

        % Equivalent ideal rate-loop commands
        q_c = (theta_c - theta) / P.tau_theta;
        p_c = (phi_c - phi) / P.tau_phi;
        r_c = 0.0;

        R.alpha(k)   = alpha;
        R.alpha_c(k) = alpha_c;
        R.alpha_n(k) = alpha_n;
        R.theta_c(k) = theta_c;
        R.phi_c(k)   = phi_c;
        R.a_n(k)     = a_n;
        R.CL(k)      = CL;
        R.CD(k)      = CD;
        R.rho(k)     = rho;
        R.q_c(k)     = q_c;
        R.p_c(k)     = p_c;
        R.r_c(k)     = r_c;
    end
end

function rho = airDensity(z, P)
% Simple exponential atmosphere
    h = max(z, 0.0);
    rho = P.rho0 * exp(-h / P.h_scale);
end

function y = saturate(u, umin, umax)
    y = min(max(u, umin), umax);
end

function y = signedLowerBound(u, boundAbs)
% Keep the sign while preventing near-zero denominators
    if abs(u) < boundAbs
        y = sign(u + eps) * boundAbs;
    else
        y = u;
    end
end

function [value, isterminal, direction] = stopEvents(~, X, P)
% Event 1: altitude reaches the ground
% Event 2: speed becomes too small
    z = X(3);
    V = X(4);

    value = [z; V - (P.V_min - 5)];
    isterminal = [1; 1];
    direction = [-1; -1];
end

function printSummary(R, P, SC)
    eV_end = P.V_ref - R.V(end);
    ean_end = P.a_n_ref - R.a_n(end);
    fprintf('Scenario: %s\n', SC.name);
    fprintf('  Final time      : %7.2f s\n', R.t(end));
    fprintf('  Final altitude  : %7.2f m\n', R.z(end));
    fprintf('  Final range x   : %7.2f m\n', R.x(end));
    fprintf('  Final speed     : %7.2f m/s   (error = %+7.2f)\n', R.V(end), eV_end);
    fprintf('  Final a_n       : %7.2f m/s^2 (error = %+7.2f)\n', R.a_n(end), ean_end);
    fprintf('  Max |alpha|     : %7.2f deg\n', max(abs(R.alpha)) * 180 / pi);
    fprintf('  Max |q_c|       : %7.2f deg/s\n\n', max(abs(R.q_c)) * 180 / pi);
end
