%% Clean workspace
clear all
clc

% Define the objective function
function cost = computeRMS(parameters)
    % Extract optimization variables
    R_val = parameters(1);      % Resistor value (Ω)
    L_val = parameters(2);      % Inductor value (H)
    C_val = parameters(3);      % Capacitor value (F)
    D_Ron = parameters(4);     % Diode forward resistance (Ω)
    Vf_Ron = parameters(5);   % Diode forward voltage (V)
    D_rev = parameters(6);     % Diode reverse resistance (Ω)
    Vf_rev = parameters(7);    % Diode reverse voltage (V)

    % Update parameters in Simulink
    set_param('rlc_full/resistor', 'R', num2str(R_val));
    set_param('rlc_full/inductor', 'L', num2str(L_val));
    set_param('rlc_full/capacitor', 'C', num2str(C_val));
    set_param('rlc_full/diode', 'Ron', num2str(D_Ron));
    set_param('rlc_full/diode', 'Vf', num2str(Vf_Ron));
    set_param('rlc_full/diode1', 'Ron', num2str(D_rev));
    set_param('rlc_full/diode1', 'Vf', num2str(Vf_rev));

    % Run the simulation
    simOut = sim('rlc_full.slx');

    % Load output signal
    % Extract Jr acceleration data from the simulation output
    load("rlc-simple/h.mat", "H2631");
    H2631 = H2631.Data; % Extract acceleration data

    % Compute RMS of the filtered acceleration
    cost = rms(H2631); % Calculate RMS acceleration (Jr metric)

    % Print tracking info
    fprintf('R: %.2f Ω, L: %.4f H, C: %.6f F, D_Ron: %.3f Ω, Vf_Ron: %.2f V, D_rev: %.3f Ω, Vf_rev: %.2f V, RMS: %.10f\n', ...
        R_val, L_val, C_val, D_Ron, Vf_Ron, D_rev, Vf_rev, cost);
end

% Define parameter bounds
lower_bounds = [1, 1e-1, 1e-2, 0.01, 0.1, 0.01, 0.1];   % [R, L, C, D_Ron, Vf_Ron, D_rev, Vf_rev]
upper_bounds = [50, 5, 1, 1, 1, 1, 1];       % [R, L, C, D_Ron, Vf_Ron, D_rev, Vf_rev]

% Initial guess
initial_guess = [8.25, 0.4125, 0.01, 0.983, 0.7, 0.983, 0.7];

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
optimal_Vf_Ron = optimal_params(5);
optimal_D_rev = optimal_params(6);
optimal_Vf_rev = optimal_params(7);

% Display results
fprintf('\n====== Optimal Parameters for rlc-diode ======\n');
fprintf('Optimal Resistance: %.2f Ω\n', optimal_R);
fprintf('Optimal Inductance: %.6f H\n', optimal_L);
fprintf('Optimal Capacitance: %.6f F\n', optimal_C);
fprintf('Optimal Diode Ron: %.3f Ω\n', optimal_D_Ron);
fprintf('Optimal Diode Forward Voltage: %.2f V\n', optimal_Vf_Ron);
fprintf('Optimal Diode Reverse Resistance: %.3f Ω\n', optimal_D_rev);
fprintf('Optimal Diode Reverse Voltage: %.2f V\n', optimal_Vf_rev);
fprintf('Minimum RMS Output: %.10f\n', optimal_rms);
fprintf('================================================\n');

% Update model one last time
set_param('rlc_full/resistor', 'R', num2str(optimal_R));
set_param('rlc_full/inductor', 'L', num2str(optimal_L));
set_param('rlc_full/capacitor', 'C', num2str(optimal_C));
set_param('rlc_full/diode', 'Ron', num2str(optimal_D_Ron));
set_param('rlc_full/diode', 'Vf', num2str(optimal_Vf_Ron));
set_param('rlc_full/diode', 'Rev', num2str(optimal_D_rev));
set_param('rlc_full/diode', 'Vf', num2str(optimal_Vf_rev));
 