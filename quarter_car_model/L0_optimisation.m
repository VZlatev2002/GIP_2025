% Define the objective function
function cost = computeRMS(parameters)
    % Update the damping coefficient in the Simulink model
    damping_value = parameters(1); % Assuming only one design parameter
    set_param('L0model/L0 (damper)', 'D', num2str(damping_value)); % Replace 'L0/damper' with your actual block path
    % Run the simulation
    simOut = sim('L0model'); % Replace with your Simulink model name
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

    % Display current values
    fprintf('Damping: %.2f Ns/m, RMS: %.6f m/s^2\n', damping_value, cost);
end

% Define initial guess and bounds for damping coefficient
lower_bound = 300;    % Minimum damping coefficient (adjust as needed)
upper_bound = 700;   % Maximum damping coefficient (adjust as needed)

% Run the optimization using fminbnd
options = optimset('Display', 'iter'); % Display optimization progress
[optimal_c, optimal_rms] = fminbnd(@computeRMS, lower_bound, upper_bound, options);

% Display results
fprintf('Optimal Damping Coefficient: %.2f Ns/m\n', optimal_c);
fprintf('Minimum RMS Acceleration: %.4f m/s^2\n', optimal_rms);

% Update the Simulink model with the optimal damping coefficient
set_param('L0model/L0 (damper)', 'D', num2str(optimal_c));
