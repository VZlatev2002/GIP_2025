%% Clean workspace
clear all
clc
warning('off','all')



%% Define the objective function
function cost = computeRMS(parameters)
    % Extract optimization variables
    Kmax = parameters(1);  % Maximum spring coefficient
    Kmin = parameters(2);  % Minimum spring coefficient
    c    = parameters(3);  % Damping coefficient
    b    = parameters(4);  % Inerter value (kg)

    % Ensure model is loaded
    model_name = 'spring_NL3';
    if ~bdIsLoaded(model_name)
        load_system(fullfile('types','spring','spring_NL3.slx'));
    end

    % Update parameters in Simulink
    set_param('spring_NL3/VDS/Kmax',        'constant', num2str(Kmax)); 
    set_param('spring_NL3/VDS/Kmin',        'constant', num2str(Kmin)); 
    set_param('spring_NL3/KL3_damper',      'D',        num2str(c)); 
    set_param('spring_NL3/KL3_inerter',     'B',        num2str(b));

    % Run the simulation
    simOut = sim(model_name);
    
    % Extract acceleration data from simulation output
    H2631 = simOut.Z2A2631;  % Access acceleration timeseries

    % Compute RMS of acceleration
    cost = rms(H2631);

    % Print tracking info
    fprintf('Kmax: %.2f, Kmin: %.2f, Damper: %.2f, Inerter: %.2f, RMS: %.6f\n', ...
            Kmax, Kmin, c, b, cost);
end

%% Define bounds for parameters
lower_bounds = [100, 10,   1,   0.1];   % [Kmax, Kmin, c, b]
upper_bounds = [20000, 20000, 500, 100];    % [Kmax, Kmin, c, b]

%% Initial guess for the parameters
initial_guess = [19000, 1900, 300, 50];    % [Kmax, Kmin, c, b]

%% Optimization options (without parallel execution)
options = optimoptions('patternsearch', ...
    'Display', 'iter', ...
    'UseCompletePoll', true, ...
    'UseCompleteSearch', true, ...
    'PollMethod', 'GPSPositiveBasis2N', ...
    'MeshTolerance', 1e-3, ...
    'UseParallel', false);   % explicitly disable parallel evaluations

%% Run optimization
[optimal_params, optimal_rms] = patternsearch(@computeRMS, initial_guess, ...
    [], [], [], [], lower_bounds, upper_bounds, [], options);

%% Extract optimized values
optimal_Kmax = optimal_params(1);
optimal_Kmin = optimal_params(2);
optimal_c    = optimal_params(3);
optimal_b    = optimal_params(4);

%% Display final optimized results
fprintf('\nOptimal Kmax: %.2f N/m\n', optimal_Kmax);
fprintf('Optimal Kmin: %.2f N/m\n', optimal_Kmin);
fprintf('Optimal Damping: %.2f Ns/m\n', optimal_c);
fprintf('Optimal Inerter: %.2f kg\n', optimal_b);
fprintf('Minimum RMS Acceleration: %.6f m/s^2\n', optimal_rms);

%% Update final Simulink model with optimized parameters
set_param('spring_NL3/VDS/Kmax',        'constant', num2str(optimal_Kmax)); 
set_param('spring_NL3/VDS/Kmin',        'constant', num2str(optimal_Kmin)); 
set_param('spring_NL3/KL3_damper',      'D',        num2str(optimal_c));
set_param('spring_NL3/KL3_inerter',     'B',        num2str(optimal_b));
