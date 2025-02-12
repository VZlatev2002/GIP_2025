clear all
function [time, acc_data] = load_acceleration_data(filename)
    % Load dataset from .mat file
    load(filename, 'Ms_acc');
    time = Ms_acc.Time;
    acc_data = Ms_acc.Data;
end

function [f, amplitude_spectrum] = compute_psd(acc_data, Fs)
    segment_length = 7500;
    overlap = segment_length * 0.8;
    nfft = 2^nextpow2(segment_length);
    window = flattopwin(segment_length);
    [psd_acc, f] = pwelch(acc_data, window, overlap, nfft, Fs);
    amplitude_spectrum = sqrt(psd_acc);
end

% Load dataset
[NL1_time, NL1_acc] = load_acceleration_data('NL5_a.mat');

% Ensure sampling frequency is consistent
dt = NL1_time(2) - NL1_time(1);
Fs = 1 / dt;

% Compute PSD and Amplitude Spectrum
[f, amplitude_NL1] = compute_psd(NL1_acc, Fs);

% Load Jr value from Simulink output
load("NL5_H.mat", "H2631");
H2631 = H2631.Data;
HRMS_NL1 = rms(H2631);
% Print Jr value
fprintf('Jr Value from Simulink Model: %.10f m/s^2\n', HRMS_NL1);