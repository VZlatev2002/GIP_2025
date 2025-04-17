% Clear the workspace
clear all
clc

% Function to perform FFT
function result = FFT_JLR(a_sprung_ratley, STEP_SIMULATION, STEP_FREQUENCY, t1, t2, f_range)
    [~, index_10s] = min(abs(a_sprung_ratley(:, 1) - t1));
    [~, index_40s] = min(abs(a_sprung_ratley(:, 1) - t2));
    x = a_sprung_ratley(index_10s:index_40s, 2);

    Fs = 1 / STEP_SIMULATION;

    w_length = 20; % Window length in seconds
    N_win = w_length / STEP_SIMULATION;
    window = flattopwin(N_win); % Flat top window
    noverlap = 0.85 * N_win;

    [SxxWelch, f] = pwelch(x, window, noverlap, f_range, Fs);
    rmsDensity = SxxWelch.^0.5;
    result = [f' rmsDensity'];
end

% Function to process data and perform FFT
function result = findFFT(file_path, STEP_SIMULATION, STEP_FREQUENCY, t1, t2, f_range)
    load(file_path);
    if exist('ans', 'var')
        time = ans.Time;
        acc_data = ans.Data;
    else
        error('Data not found in the file.');
    end

    tsimnew = [time(1):STEP_SIMULATION:time(end)]'; % Ensure constant time step
    msanew = interp1(time, acc_data, tsimnew);

    your_acceleration_data = [tsimnew, msanew];

    % Process FFT
    FFT_result = FFT_JLR(your_acceleration_data, STEP_SIMULATION, STEP_FREQUENCY, t1, t2, f_range);
    result = FFT_result;
end

% Simulation Parameters
STEP_SIMULATION = 0.0001; % Time step
STEP_FREQUENCY = 0.1; % Frequency resolution
t1 = 0; % Start time for FFT
t2 = 100; % End time for FFT
f_range = [0:STEP_FREQUENCY:20]; % Frequency range for FFT

% Load and process data for each "a" experiment
result_rlc_series = findFFT('data_elec/a_rlc_s.mat', STEP_SIMULATION, STEP_FREQUENCY, t1, t2, f_range);
result_rlc_parallel = findFFT('data_elec/a_rlc_p.mat', STEP_SIMULATION, STEP_FREQUENCY, t1, t2, f_range);
result_a1_parallel = findFFT('data_elec/a_1_parallel.mat', STEP_SIMULATION, STEP_FREQUENCY, t1, t2, f_range);
result_a1_series = findFFT('data_elec/a_1_series.mat', STEP_SIMULATION, STEP_FREQUENCY, t1, t2, f_range);
result_a2_parallel = findFFT('data_elec/a_2_parallel.mat', STEP_SIMULATION, STEP_FREQUENCY, t1, t2, f_range);
result_a2_series = findFFT('data_elec/a_2_series.mat', STEP_SIMULATION, STEP_FREQUENCY, t1, t2, f_range);

% Plot the results
figure;
hold on;
plot(result_rlc_series(:, 1), result_rlc_series(:, 2), 'LineWidth', 1.5, 'DisplayName', 'RLC-S(default)');
plot(result_rlc_parallel(:, 1), result_rlc_parallel(:, 2), 'LineWidth', 1.5, 'DisplayName', 'RLC-P');
plot(result_a1_parallel(:, 1), result_a1_parallel(:, 2), 'LineWidth', 1.5, 'DisplayName', '1DP');
plot(result_a1_series(:, 1), result_a1_series(:, 2), 'LineWidth', 1.5, 'DisplayName', '1DS');
plot(result_a2_parallel(:, 1), result_a2_parallel(:, 2), 'LineWidth', 1.5, 'DisplayName', '2DP');
plot(result_a2_series(:, 1), result_a2_series(:, 2), 'LineWidth', 1.5, 'DisplayName', '2DS');

% Add legend
legend show;

% Add labels and title
xlabel('Frequency (Hz)');
ylabel('Sprung mass acceleration (m/s^2)');
title('Power Spectral density of Acceleration');

% Hold off to finalize the plot
hold off;
