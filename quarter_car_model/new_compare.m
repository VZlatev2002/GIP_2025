%% Main Script to Compare L0 and L1 Acceleration Spectra
clc; clear; close all;


function [time, acc_data] = load_acceleration_data(filename)
    % Load dataset from .mat file
    load(filename, 'Ms_acc');
    time = Ms_acc.Time;       % Extract time vector
    acc_data = Ms_acc.Data;   % Extract acceleration data
end

function [f, amplitude_spectrum] = compute_psd(acc_data, Fs)
    % Define pwelch parameters
    segment_length = 5000; % 5s segment length (5,000 samples per segment)
    overlap = segment_length * 0.5; % 50% overlap (2,500 samples)
    nfft = 2^nextpow2(segment_length); % Power of 2 for FFT efficiency
    window = flattopwin(segment_length); % Flat Top Window for accurate peaks

    % Compute Power Spectral Density (PSD)
    [psd_acc, f] = pwelch(acc_data, window, overlap, nfft, Fs);

    % Convert PSD to Amplitude Spectrum
    amplitude_spectrum = sqrt(psd_acc);

end

function plot_acceleration_spectrum(f, amplitude_L0, amplitude_L1,amplitude_NL0,amplitude_NL1)
    figure;
    plot(f, amplitude_L0, 'b', 'LineWidth', 1.5, 'DisplayName', 'L0 Model');
    hold on;
    plot(f, amplitude_L1, 'r', 'LineWidth', 1.5, 'DisplayName', 'L1 Model');
    plot(f, amplitude_NL0, 'black', 'LineWidth', 1.5, 'DisplayName', 'NL0 Model');
    plot(f, amplitude_NL1, 'green', 'LineWidth', 1.5, 'DisplayName', 'NL1 Model');
    xlabel('Frequency (Hz)');
    ylabel('Acceleration Amplitude (m/s²)');
    title('Acceleration vs Frequency Spectrum (L0 vs L1) - Flat Top Window');
    grid on;
    xlim([0.1 20]);  % Focus on relevant frequency range
    legend;
    hold off;
end

function [J1] = cost(time,acc_data)
    % Define the transfer function
    numerator = [50, 500];  % Transfer function numerator
    denominator = [1, 50, 1200]; % Transfer function denominator

    % Apply the transfer function to the acceleration data
    filtered_acceleration = lsim(tf(numerator, denominator), acc_data, time);

    % Compute RMS of the filtered acceleration
    J1 = sqrt(mean(filtered_acceleration.^2)); % Calculate RMS acceleration
end

% Function to find two peak values: one below 4Hz and one above 4Hz
function [peak_low, peak_high] = find_two_peaks(f, amplitude)
    % Split frequency range into below 4Hz and above 4Hz
    idx_low = f < 4;
    idx_high = f >= 4;
    
    % Find the maximum peak in each region
    [peak_low, idx1] = max(amplitude(idx_low));  
    [peak_high, idx2] = max(amplitude(idx_high));

    % Get corresponding frequencies
    freq_low = f(idx_low);
    freq_high = f(idx_high);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load datasets
[L0_time, L0_acc] = load_acceleration_data('L0acc_data100s.mat');
[L1_time, L1_acc] = load_acceleration_data('L1acc_data100s.mat');
[NL0_time, NL0_acc] = load_acceleration_data('NL0acc_data100s.mat');
[NL1_time, NL1_acc] = load_acceleration_data('NL1acc_data100s.mat');

% Ensure sampling frequency is consistent
dt = L0_time(2) - L0_time(1);  % Time step
Fs = 1 / dt;                   % Sampling frequency (1000 Hz)

% Compute PSD and Amplitude Spectrum for all models
[f, amplitude_L0] = compute_psd(L0_acc, Fs);
[~, amplitude_L1] = compute_psd(L1_acc, Fs); % Frequency vector `f` is the same
[~, amplitude_NL0] = compute_psd(NL0_acc, Fs); % NL0 model
[~, amplitude_NL1] = compute_psd(NL1_acc, Fs); % NL0 model

% Compute RMS for all models
J1_L0 = cost(L0_time, L0_acc);
J1_L1 = cost(L1_time, L1_acc);
J1_NL0 = cost(NL0_time, NL0_acc);
J1_NL1 = cost(NL1_time, NL1_acc);

% Compute percentage improvement in RMS
J1_L1_improvement = ((J1_L0 - J1_L1) / J1_L0) * 100;
J1_NL0_improvement = ((J1_L0 - J1_NL0) / J1_L0) * 100;
J1_NL1_improvement = ((J1_L0 - J1_NL1) / J1_L0) * 100;

% Print RMS values
fprintf('J1 (L0 Model): %.6f m/s^2\n', J1_L0);
fprintf('J1 (L1 Model): %.6f m/s^2\n', J1_L1);
fprintf('J1 (NL0 Model): %.6f m/s^2\n', J1_NL0);
fprintf('J1 (NL1 Model): %.6f m/s^2\n', J1_NL1);
fprintf('Improvement in J1 (L1 vs L0): %.2f%%\n', J1_L1_improvement);
fprintf('Improvement in J1 (NL0 vs L0): %.2f%%\n', J1_NL0_improvement);
fprintf('Improvement in J1 (NL1 vs L0): %.2f%%\n', J1_NL1_improvement);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Acceleration vs Frequency for all models
plot_acceleration_spectrum(f, amplitude_L0, amplitude_L1, amplitude_NL0,amplitude_NL1);

% Find two peak values for each model
[peak_L0_low, peak_L0_high] = find_two_peaks(f, amplitude_L0);
[peak_L1_low, peak_L1_high] = find_two_peaks(f, amplitude_L1);
[peak_NL0_low, peak_NL0_high] = find_two_peaks(f, amplitude_NL0);
[peak_NL1_low, peak_NL1_high] = find_two_peaks(f, amplitude_NL1);

% Compute percentage improvement in peaks relative to L0
peak_L1_improvement_low = ((peak_L0_low - peak_L1_low) / peak_L0_low) * 100;
peak_L1_improvement_high = ((peak_L0_high - peak_L1_high) / peak_L0_high) * 100;
peak_NL0_improvement_low = ((peak_L0_low - peak_NL0_low) / peak_L0_low) * 100;
peak_NL0_improvement_high = ((peak_L0_high - peak_NL0_high) / peak_L0_high) * 100;
peak_NL1_improvement_low = ((peak_L0_low - peak_NL1_low) / peak_L0_low) * 100;
peak_NL1_improvement_high = ((peak_L0_high - peak_NL1_high) / peak_L0_high) * 100;

% Print peak values and improvements
fprintf('\nPeak Analysis:\n');
fprintf('L0 Peaks: Low: %.6f, High: %.6f\n', peak_L0_low, peak_L0_high);
fprintf('L1 Peaks: Low: %.6f (%.2f%%), High: %.6f (%.2f%%)\n', ...
    peak_L1_low, peak_L1_improvement_low, peak_L1_high, peak_L1_improvement_high);
fprintf('NL0 Peaks: Low: %.6f (%.2f%%), High: %.6f (%.2f%%)\n', ...
    peak_NL0_low, peak_NL0_improvement_low, peak_NL0_high, peak_NL0_improvement_high);
fprintf('NL1 Peaks: Low: %.6f (%.2f%%), High: %.6f (%.2f%%)\n', ...
    peak_NL1_low, peak_NL1_improvement_low, peak_NL1_high, peak_NL1_improvement_high);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load("NL1_middle_v.mat",'middle_v')
load("NL1_Ms_v.mat",'Ms_v')

time = Ms_v.Time;       % Extract time vector
Ms_v = Ms_v.Data;   % Extract acceleration data
middle_v = middle_v.Data;

relative_VDD_v = middle_v - Ms_v;
max_relative_v = max(relative_VDD_v);
Cmax = 740;
max_force = max_relative_v * Cmax;

fprintf('Max force to VDD: %.2fPa\n', max_force);



