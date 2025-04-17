%% Clean workspace
clear all
clc

% Define the objective function
function cost = computeRMS(parameters)
    % Extract optimization variables
    R_A = parameters(1);  % Low resistance state (ON)
    R_B = parameters(2);  % High resistance state (OFF)

    % Update parameters in Simulink (Memristor resistance values)
    set_param('new_electrical/custom_memristor', 'Roff', num2str(R_A)); 
    set_param('new_electrical/custom_memristor', 'Ron', num2str(R_B));

    % Run the simulation
    simOut = sim('new_electrical.slx');

    % Extract Jr acceleration data from the simulation output
    load("electrical_simplified/h.mat", "H2631");
    H2631 = H2631.Data; % Extract acceleration data

    % Compute RMS of the filtered acceleration
    cost = rms(H2631); % Calculate RMS acceleration (Jr metric)

    % Print current values for tracking
    fprintf('R_A: %.2f Ω, R_B: %.2f Ω, RMS: %.10f\n', R_A, R_B, cost);
end

% Define bounds for the parameters
lower_bounds = [10, 10];      % [Min resistance (Ω), Max resistance (Ω)]
upper_bounds = [100, 20];   % [Min resistance (Ω), Max resistance (Ω)]

% Initial guess for the parameters
initial_guess = [50, 15];  % [R_A (low), R_B (high)]

% Set optimization options
options = optimoptions('patternsearch', 'Display', 'iter', ...
        'UseCompletePoll', true, 'UseCompleteSearch', true, ...
        'PollMethod', 'GPSPositiveBasis2N', 'MeshTolerance', 1e-3);

% Run the optimization using pattern search
[optimal_params, optimal_rms] = patternsearch(@computeRMS, initial_guess, [], [], [], [], ...
                                              lower_bounds, upper_bounds, [], options);

% Extract optimized values
optimal_R_A = optimal_params(1);
optimal_R_B = optimal_params(2);
    
% Display results
fprintf('\n====== Optimal Parameters for new_electrical ======\n');
fprintf('Optimal R_A: %.2f Ω\n', optimal_R_A);
fprintf('Optimal R_B: %.2f Ω\n', optimal_R_B);
fprintf('Minimum RMS Acceleration: %.10f m/s^2\n', optimal_rms);
fprintf('===================================================\n');

% Update the Simulink model with the optimized values
set_param('new_electrical/custom_memristor', 'Roff', num2str(optimal_R_A));
set_param('new_electrical/custom_memristor', 'Ron', num2str(optimal_R_B));
