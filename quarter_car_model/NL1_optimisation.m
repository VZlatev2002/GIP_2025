% Define the objective function
function cost = computeRMS(parameters)
    % Extract optimization variables
    Cmax = parameters(1);  % Maximum damping coefficient
    Cmin = parameters(2);  % Minimum damping coefficient
    ks = parameters(3);    % Spring stiffness (N/m)
    b = parameters(4);     % Inerter value (kg)

    % Update parameters in Simulink
    set_param('NL1model/NL1 Cmax', 'Gain', num2str(Cmax)); 
    set_param('NL1model/NL1 Cmin', 'Gain', num2str(Cmin)); 
    set_param('NL1model/NL1 (spring)', 'spr_rate', num2str(ks)); 
    set_param('NL1model/NL1 (inerter)', 'B', num2str(b));

    % Run the simulation
    simOut = sim('NL1model');

    % Extract acceleration data from the simulation output
    acc_data = simOut.acc.Data; % Access the acceleration data (timeseries object)
    time1 = simOut.acc.Time;    % Extract corresponding time vector

    % Define the transfer function
    numerator = [50, 500];  % Transfer function numerator
    denominator = [1, 50, 1200]; % Transfer function denominator

    % Apply the transfer function to the acceleration data
    filtered_acceleration = lsim(tf(numerator, denominator), acc_data, time1);

    % Compute RMS of the filtered acceleration
    cost = sqrt(mean(filtered_acceleration.^2)); % Calculate RMS acceleration

    % Print current values for tracking
    fprintf('Cmax: %.2f, Cmin: %.2f, Spring: %.2f, Inerter: %.2f, RMS: %.6f\n', ...
            Cmax, Cmin, ks, b, cost);
end

% Define bounds for the parameters
lower_bounds = [100, 1, 10000, 0.1];   % [Cmax, Cmin, Spring stiffness, Inerter]
upper_bounds = [4000, 100, 100000, 50]; % [Cmax, Cmin, Spring stiffness, Inerter]

% Initial guess for the parameters
initial_guess = [500, 50, 20000, 5]; % [Cmax, Cmin, Spring stiffness, Inerter]

% Set optimization options
options = optimoptions('patternsearch', 'Display', 'iter', ...
    'UseCompletePoll', true, 'UseCompleteSearch', true, ...
    'PollMethod', 'GPSPositiveBasis2N', 'MeshTolerance', 1e-3);

% Run the optimization using pattern search
[optimal_params, optimal_rms] = patternsearch(@computeRMS, initial_guess, [], [], [], [], ...
                                              lower_bounds, upper_bounds, [], options);

% Extract optimized values
optimal_Cmax = optimal_params(1);
optimal_Cmin = optimal_params(2);
optimal_ks = optimal_params(3);
optimal_b = optimal_params(4);

% Display results
fprintf('\nOptimal Cmax: %.2f Ns/m\n', optimal_Cmax);
fprintf('Optimal Cmin: %.2f Ns/m\n', optimal_Cmin);
fprintf('Optimal Spring Stiffness: %.2f N/m\n', optimal_ks);
fprintf('Optimal Inerter: %.2f kg\n', optimal_b);
fprintf('Minimum RMS Acceleration: %.6f m/s^2\n', optimal_rms);

% Update the Simulink model with the optimized values
set_param('NL1model/NL1 Cmax', 'Gain', num2str(optimal_Cmax));
set_param('NL1model/NL1 Cmin', 'Gain', num2str(optimal_Cmin));
set_param('NL1model/NL1 (spring)', 'spr_rate', num2str(optimal_ks));
set_param('NL1model/NL1 (inerter)', 'B', num2str(optimal_b));
