clc;
clear;
close all;

% MULTI-ATTITUDE TRACKING MANEUVERS

% Define initial conditions and parameters based on the paper
omega0 = [0.001; -0.005; 0.001];   % Initial angular velocity
rho0 = [1; 1; -1];                  % Initial attitude components
rho_d = @(t) [sin(pi*t/50); -sin(pi*t/50); 0.5*cos(pi*t/50)]; % Reference vector

% Parameters for J0 and Delta J (variation in inertia)
J0 = diag([87.212, 86.067, 114.562]);  % Nominal inertia matrix
deltaJ_max = diag([8.7212, 4.3034, 17.1843]); % Maximum allowable variation
J = J0 + deltaJ_max;  % Perturbed inertia matrix

% Control parameters
lambda = 0.5; % Convergence rate of sliding mode control
epsilon = 0.05; % Saturation function parameter for SMC

% Simulation time
T = 200; % Total time in seconds (adjust as needed)
dt = 0.01; % Time step in seconds (adjust as needed)
time = 0:dt:T; % Time vector
N = length(time); % Number of time steps

d = zeros(3, N);
d(1, :) = -0.005 * sin(time);
d(2, :) = 0.005 * sin(time);
d(3, :) = -0.005 * sin(time);


% Initialize state variables for storing results
rho = zeros(3, length(time));
omega = zeros(3, length(time));
tau_control = zeros(3, length(time));

rho_d_values = zeros(3, length(time));
omega_cap = zeros(3, length(time));
T_rho_omega = zeros(3, length(time));

% Set initial conditions
rho(:,1) = rho0;
omega(:,1) = omega0;
rho_d_dot0 = [pi/50; -pi/50; 0];
% omega_cap(:,1) = 0.5 * (eye(3) - skew_symmetric(rho0)) * rho_d_dot0;
omega_cap(:,1) = (2/(1 + rho0'*rho0)) * (eye(3) - skew_symmetric(rho0)) * rho_d_dot0;

for i = 1:length(time)-1
    t = time(i);
    
    % Desired reference and its derivative
    rho_d_t = rho_d(t);
    rho_d_dot = (rho_d(t + dt) - rho_d(t)) / dt; % Approximate derivative
    rho_d_values(:, i) = rho_d(t);
    
    % Compute the transformation matrix T(rho)
    rho_t = rho(:,i);
    rho_cross = skew_symmetric(rho_t);
    T_rho = 0.5 * (eye(3) + (rho_t * rho_t') + rho_cross);
    rho_dot = T_rho * omega(:, i);
    T_rho_inv = 2 / (1 + rho_t' * rho_t) * (eye(3) - rho_cross);

    % Sliding vector s
    s = omega(:,i) - omega_cap(:,i) + lambda * (rho_t - rho_d_t);

    rho_d_2dot = (rho_d(t + 2*dt) - 2*rho_d(t + dt) + rho_d(t)) / (dt^2);

    rho_T_rho = rho_t' * rho_t;                     % Scalar: rho^T * rho
    T_rho_omega(:,i) = T_rho * omega(:,i);                % Vector: T(rho) * omega
    
    term1 = -2 / (1 + rho_T_rho) * skew_symmetric(T_rho_omega(:,i));
   
    rho_T_rho_omega = rho_t' * T_rho_omega(:,i);

    term2 = -4 * rho_T_rho_omega / (1 + rho_T_rho)^2 * (eye(3) - skew_symmetric(rho_t));
    
    dT_inv_dt = term1 + term2;

    omega_cap_dot = dT_inv_dt * rho_d_dot + T_rho_inv * rho_d_2dot;

    k1 = abs((deltaJ_max(2,2) - deltaJ_max(3,3)) * omega(2,i) * omega(3,i)) + deltaJ_max(1,1) * abs(omega_cap_dot(1)) ...
        + lambda * deltaJ_max(1,1) * abs(T_rho_omega(1,i) - rho_d_dot(1)) + abs(0.005) + 1;

    k2 = abs((deltaJ_max(3,3) - deltaJ_max(1,1)) * omega(3,i) * omega(1,i)) + deltaJ_max(2,2) * abs(omega_cap_dot(2)) ...
        + lambda * deltaJ_max(2,2) * abs(T_rho_omega(2,i) - rho_d_dot(2)) + abs(0.005) + 1;

    k3 = abs((deltaJ_max(1,1) - deltaJ_max(2,2)) * omega(1,i) * omega(2,i)) + deltaJ_max(3,3) * abs(omega_cap_dot(3)) ...
        + lambda * deltaJ_max(3,3) * abs(T_rho_omega(3,i) - rho_d_dot(3)) + abs(0.005) + 1;
     
    h0 = J0 * omega(:, i);
        
    % control torque computation
    tau = - skew_symmetric(h0) * omega(:,i) + J0 * omega_cap_dot ...
          - J0 * lambda * (T_rho * omega(:,i) - rho_d_dot) ...
          - K_sign(s, epsilon, k1, k2, k3);
    tau_control(:,i) = tau;

    % current disturbance
    d_current = d(:, i);                     % 3x1 vector for the current time step
   
    H_omega_term = skew_symmetric(h0) * omega(:, i);  % H(omega)*omega
    tau_plus_d = tau + d_current;             % Control + disturbance

    omega_dot = J \ (tau_plus_d + H_omega_term); 

    omega(:, i+1) = omega(:, i) + omega_dot * dt;

    rho(:,i+1) = rho(:,i) + rho_dot * dt;

    omega_cap(:,i+1) = omega_cap(:,i) + omega_cap_dot * dt;
end


figure(1);
hold on;
plot(time, rho(1,:), 'r', time, rho(2,:), 'g', time, rho(3,:), 'b', 'LineWidth',1.5);
plot(time, rho_d_values(1,:), 'r--', time, rho_d_values(2,:), 'g--', time, rho_d_values(3,:), 'b--', 'LineWidth',1.5);
hold off;
grid on;
ylim([-2, 1.5]);
xlabel('Time (s)');
ylabel('\rho components');
title('Attitude Components \rho_i, i = 1, 2, 3');
legend('\rho_1', '\rho_2', '\rho_3','Location','best');

figure(2);
hold on;
plot(time, omega(1,:), 'r', time, omega(2,:), 'g', time, omega(3,:), 'b','LineWidth',1.5);
plot(time, omega_cap(1,:), 'r--', time, omega_cap(2,:), 'g--', time, omega_cap(3,:), 'b--','LineWidth',1.5);
hold off;
grid on;
xlabel('Time (s)');
ylabel('\omega components');
title('Angular Velocities \omega_i, i = 1, 2, 3');
legend('\omega_1', '\omega_2', '\omega_3','Location','best');

figure(3);
plot(time, tau_control(1,:), 'r', time, tau_control(2,:), 'g', time, tau_control(3,:), 'b','LineWidth',1.5);
grid on;
xlabel('Time (s)');
ylabel('\tau components');
title('Control Torques \tau_i, i = 1, 2, 3');
legend('\tau_1', '\tau_2', '\tau_3','Location','best');


function K_s = K_sign(s, epsilon, k1, k2, k3)
    K_s = [k1; k2; k3] .* sat(s, epsilon);
end

function sat_val = sat(s, epsilon)
    sat_val = min(1, abs(s / epsilon)) .* sign(s);
end


function S = skew_symmetric(v)
    S = [   0    -v(3)   v(2);
           v(3)    0    -v(1);
          -v(2)   v(1)    0  ];
end
