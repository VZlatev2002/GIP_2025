%% === Electrical Circuit Optimization Script ===
% This script optimizes three resistors (R1, R2, R3), a diode (D1),
% and a capacitor (C1) in a Simulink model to minimize RMS voltage ripple.

clc; clear; close all;

% Define parameter bounds
lower_bounds = [1, 0.0001, 1, 0.2, 1e-9];      % Min values: [R1 (Ω), R2 (Ω), R3 (Ω), D1 Voltage (V), C1 (F)]
upper_bounds = [1000, 1000, 1000, 1.0, 1e-3]; % Max values: [R1 (Ω), R2 (Ω), R3 (Ω), D1 Voltage (V), C1 (F)]

% Initial guess for parameters
initial_guess = [10, 0.001, 10, 0.3, 1e-5]; % Example: R1=10Ω, R2=0.001Ω, R3=10Ω, Diode=0.3V, Cap=1μF

% Set optimization options
options = optimoptions('patternsearch', 'Display', 'iter', ...
    'UseCompletePoll', true, 'UseCompleteSearch', true, ...
    'PollMethod', 'GPSPositiveBasis2N', 'MeshTolerance', 1e-3);

% Run the optimization using pattern search
[optimal_params, optimal_rms] = patternsearch(@computeRMS, initial_guess, [], [], [], [], ...
                                              lower_bounds, upper_bounds, [], options);

% Extract optimized values
optimal_R1 = optimal_params(1);
optimal_R2 = optimal_params(2);
optimal_R3 = optimal_params(3);
optimal_D1 = optimal_params(4);
optimal_C1 = optimal_params(5);

% Display results
fprintf('\n=== Optimization Results ===\n');
fprintf('Optimal R1: %.4f Ω\n', optimal_R1);
fprintf('Optimal R2: %.4f Ω\n', optimal_R2);
fprintf('Optimal R3: %.4f Ω\n', optimal_R3);
fprintf('Optimal Diode Voltage: %.3f V\n', optimal_D1);
fprintf('Optimal Capacitor C1: %.6e F\n', optimal_C1);
fprintf('Minimum RMS Ripple: %.10f V\n', optimal_rms);

% Update Simulink model with optimized values
set_param('VDD_NL0_electrical/R1', 'R1', num2str(optimal_R1));
set_param('VDD_NL0_electrical/R2', 'R2', num2str(optimal_R2));
set_param('VDD_NL0_electrical/R3', 'R3', num2str(optimal_R3));
set_param('VDD_NL0_electrical/D1', 'D1', num2str(optimal_D1));
set_param('VDD_NL0_electrical/C1', 'C1', num2str(optimal_C1));

%% === Objective Function: Compute RMS Ripple ===
function cost = computeRMS(parameters)
    % Extract optimization variables
    R1_val = parameters(1); % Resistor R1 (Ω)
    R2_val = parameters(2); % Resistor R2 (Ω)
    R3_val = parameters(3); % Resistor R3 (Ω)
    D1_voltage = parameters(4); % Diode Forward Voltage (V)
    C1_val = parameters(5); % Capacitor C1 (F)

    % Update parameters in the Simulink model
    set_param('VDD_NL0_electrical/R1', 'R', num2str(R1_val));
    set_param('VDD_NL0_electrical/R2', 'R', num2str(R2_val));
    set_param('VDD_NL0_electrical/R3', 'R', num2str(R3_val));
    set_param('VDD_NL0_electrical/D1', 'Vf', num2str(D1_voltage));
    set_param('VDD_NL0_electrical/C1', 'c', num2str(C1_val));

    % Run the Simulink simulation
    simOut = sim('VDD_NL0_electrical'); % Ensure this is your correct Simulink model name

    % Extract simulation output variable
    A = simOut.Z2A2631; % Your variable from the Simulink model

    % Compute RMS using your method
    cost = rms(A);

    % Print current iteration values
    fprintf('R1: %.4f Ω, R2: %.4f Ω, R3: %.4f Ω, D1 Voltage: %.3f V, C1: %.6e F, Jr: %.10f V\n', ...
        R1_val, R2_val, R3_val, D1_voltage, C1_val, cost);
end
