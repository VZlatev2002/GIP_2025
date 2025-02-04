% Define the objective function
function cost = computeRMS(parameters)
    % Extract optimization variables
    Cmax = parameters(1);  % Maximum damping coefficient
    Cmin = parameters(2);  % Minimum damping coefficient

    % Update the damping coefficients in Simulink
    set_param('NL0model/NL0 Cmax', 'Gain', num2str(Cmax)); 
    set_param('NL0model/NL0 Cmin', 'Gain', num2str(Cmin)); 

    % Run the simulation
    simOut = sim('NL0model'); % Replace 'L1model' with your Simulink model name

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

    % Print current values
    fprintf('Cmax: %.2f, Cmin: %.2f, RMS: %.6f\n', Cmax, Cmin, cost);
end
% Define bounds for the parameters
lower_bounds = [10, 1];    % Lower bounds for [Cmax, Cmin]
upper_bounds = [2000, 1000]; % Upper bounds for [Cmax, Cmin]

% Initial guess for the parameters
initial_guess = [500, 50]; % Initial values for [Cmax, Cmin]

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

% Display results
fprintf('Optimal Cmax: %.2f Ns/m\n', optimal_Cmax);
fprintf('Optimal Cmin: %.2f Ns/m\n', optimal_Cmin);
fprintf('Minimum RMS Acceleration: %.6f m/s^2\n', optimal_rms);

% Update the Simulink model with the optimized values
set_param('NL0model/NL0 Cmax', 'Gain', num2str(optimal_Cmax));
set_param('NL0model/NL0 Cmin', 'Gain', num2str(optimal_Cmin));
