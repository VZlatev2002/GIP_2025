% Define the objective function
function cost = computeRMS(parameters)
    % Update the damping coefficient in the Simulink model
    damping_value = parameters(1); % Assuming only one design parameter
    spring_value = parameters(2);   % Spring stiffness
    inertance_value = parameters(3); % Inertance

    set_param('L1model/L1 (damper)', 'D', num2str(damping_value)); % Replace 'L0/damper' with your actual block path
    set_param('L1model/L1 (spring)', 'spr_rate', num2str(spring_value));  % Update spring stiffness
    set_param('L1model/L1 (inerter)', 'B', num2str(inertance_value)); % Update inertance
    % Run the simulation
    simOut = sim('L1model'); % Replace with your Simulink model name
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

    % Print the parameters and RMS
    fprintf('Damping: %.2f, Stiffness: %.2f, Inertance: %.2f, RMS: %.6f\n', ...
            damping_value, spring_value, inertance_value, cost);
    
end


% Define bounds for the parameters
lower_bounds = [100, 1000, 0.001];  % Lower bounds for [damping, stiffness, inertance]
upper_bounds = [2000, 20000, 50]; % Upper bounds for [damping, stiffness, inertance]

% Initial guess for the parameters
initial_guess = [100, 12000, 1]; % [damping, stiffness, inertance]

% Run the optimization using fmincon
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
    'OptimalityTolerance', 1e-3, 'ConstraintTolerance', 1e-3);

% Run the global optimization
[optimal_params, optimal_rms] = patternsearch(@computeRMS, initial_guess, [], [], [], [], ...
                                              lower_bounds, upper_bounds, [], options);

% Display results
fprintf('Optimal Damping Coefficient: %.2f Ns/m\n', optimal_params(1));
fprintf('Optimal Spring Stiffness: %.2f N/m\n', optimal_params(2));
fprintf('Optimal Inertance: %.2f kg\n', optimal_params(3));
fprintf('Minimum RMS Acceleration: %.6f m/s^2\n', optimal_rms);

% Update the Simulink model with the optimal parameters
set_param('L1model/L1 (damper)', 'D', num2str(optimal_params(1)));
set_param('L1model/L1 (spring)', 'spr_rate', num2str(optimal_params(2)));
set_param('L1model/L1 (inerter)', 'B', num2str(optimal_params(3)));




