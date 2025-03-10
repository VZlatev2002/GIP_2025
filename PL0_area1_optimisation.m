% Define the objective function
function cost = computeRMS(area1)
    % Update the restriction area in the Simulink model
    %set_param('VDD_NL0/Local Restriction (IL)', 'restriction_area', num2str(area1));
    set_param('VDD_NL0/Orifice (IL)', 'orifice_area_constant', num2str(area1));
    % Run the simulation
    simOut = sim('VDD_NL0.slx'); % Replace with your Simulink model name

    % Extract acceleration data from the simulation output
    A = simOut.Z2A2631; % Access the acceleration data (timeseries object)

    % Compute RMS of the filtered acceleration
    cost = rms(A); % Calculate RMS acceleration

    % Display current values
    fprintf('Area1: %.10f Ns/m, RMS: %.6f m/s^2\n', area1, cost);
end

% Define bounds for the parameter
lower_bound = 1e-5;    % Lower bound for Area1
upper_bound = 20e-5;  % Upper bound for Area1

% Initial guess for the parameter
initial_guess = 2e-5; % Initial value for Area1

% Set optimization options
options = optimoptions('patternsearch', 'Display', 'iter', ...
    'UseCompletePoll', true, 'UseCompleteSearch', true, ...
    'PollMethod', 'GPSPositiveBasis2N', 'MeshTolerance', 1e-3);

% Run the optimization using pattern search
[optimal_area, optimal_rms] = patternsearch(@computeRMS, initial_guess, [], [], [], [], ...
                                            lower_bound, upper_bound, [], options);

% Display results
fprintf('Optimal Area1: %.10f Ns/m\n', optimal_area);
fprintf('Minimum RMS Acceleration: %.6f m/s^2\n', optimal_rms);

% Update the Simulink model with the optimized value
%set_param('VDD_NL0/Local Restriction (IL)', 'restriction_area', num2str(optimal_area));
set_param('VDD_NL0/Orifice (IL)', 'orifice_area_constant', num2str(optimal_area));
%get_param('VDD_NL0/Orifice (IL)','ObjectParameters')