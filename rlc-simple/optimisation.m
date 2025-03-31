%% Clean workspace
clear all
clc

% Define the objective function
function cost = computeRMS(parameters)
    % Extract optimization variables
    R_val = parameters(1);      % Resistor value (Ω)
    L_val = parameters(2);      % Inductor value (H)
    C_val = parameters(3);      % Capacitor value (F)
    D_Ron = parameters(4);      % Diode forward resistance (Ω)

    % Update parameters in Simulink
    set_param('rlc_full/resistor', 'R', num2str(R_val));
    set_param('rlc_full/inductor', 'l', num2str(L_val));
    set_param('rlc_full/capacitor', 'c', num2str(C_val));
    set_param('rlc_full/diode', 'Ron', num2str(D_Ron));

    % Run the simulation
    simOut = sim('rlc_full.slx');

    % Load output signal
    % Extract Jr acceleration data from the simulation output
    load("rlc-simple/h.mat", "H2631");
    H2631 = H2631.Data; % Extract acceleration data

    % Compute RMS of the filtered acceleration
    cost = rms(H2631); % Calculate RMS acceleration (Jr metric)

    % Print tracking info
    fprintf('R: %.2f Ω, L: %.4f H, C: %.6f F, D_Ron: %.3f Ω, RMS: %.10f\n', ...
        R_val, L_val, C_val, D_Ron, cost);
end

% Define parameter bounds
lower_bounds = [1, 1e-1, 1e-2, 0.01];   % [R, L, C, D_Ron]
upper_bounds = [50, 5, 1, 1];       % [R, L, C, D_Ron]

% Initial guess
initial_guess = [8.25, 0.4125, 0.01, 0.983];

% Optimization options
options = optimoptions('patternsearch', 'Display', 'iter', ...
    'UseCompletePoll', true, 'UseCompleteSearch', true, ...
    'PollMethod', 'GPSPositiveBasis2N', 'MeshTolerance', 1e-3);

% Run optimization
[optimal_params, optimal_rms] = patternsearch(@computeRMS, initial_guess, [], [], [], [], ...
                                              lower_bounds, upper_bounds, [], options);

% Extract optimized values
optimal_R = optimal_params(1);
optimal_L = optimal_params(2);
optimal_C = optimal_params(3);
optimal_D_Ron = optimal_params(4);

% Display results
fprintf('\n====== Optimal Parameters for rlc-diode ======\n');
fprintf('Optimal Resistance: %.2f Ω\n', optimal_R);
fprintf('Optimal Inductance: %.6f H\n', optimal_L);
fprintf('Optimal Capacitance: %.6f F\n', optimal_C);
fprintf('Optimal Diode Ron: %.3f Ω\n', optimal_D_Ron);
fprintf('Minimum RMS Output: %.10f\n', optimal_rms);
fprintf('================================================\n');

% Update model one last time
set_param('rlc_full/resistor', 'resistor', num2str(optimal_R));
set_param('rlc_full/inductor', 'inductor', num2str(optimal_L));
set_param('rlc_full/capacitor', 'capacitor', num2str(optimal_C));
set_param('rlc_full/diode', 'Ron', num2str(optimal_D_Ron));
