% Parameters
V = 25; % Vehicle speed (m/s)
kappa = 5e-7; % Roughness parameter (m^3/cycle)
f_min = 0.1; % Minimum frequency (Hz)
f_max = 100; % Maximum frequency (Hz)
N = 2^14; % Number of samples
fs = 1000; % Sampling frequency (Hz)
df = fs / N; % Frequency resolution
time = (0:N-1) / fs; % Time vector

% Frequency vector
freq = (0:N/2) * df; % Frequency vector (Hz), size = N/2+1

% Road velocity PSD
S_vel = @(f) 2 * pi * kappa * V; % Velocity PSD is constant
psd_values = S_vel(freq); % Evaluate PSD for all frequencies

% Add Gaussian white noise
white_noise = randn(1, N); % Generate Gaussian white noise
noise_spectrum = fft(white_noise); % Convert to frequency domain, size = N

% Shape the spectrum with PSD
amplitude = sqrt([psd_values, fliplr(psd_values(2:end-1))] * df); % Mirror PSD to match size N
random_phases = exp(1j * 2 * pi * rand(1, N/2 + 1)); % Random phases for positive frequencies
random_phases_full = [random_phases, conj(fliplr(random_phases(2:end-1)))]; % Symmetric phases, size N
shaped_spectrum = amplitude .* random_phases_full; % Create shaped spectrum, size N

% Add noise to the shaped spectrum
final_spectrum = noise_spectrum .* shaped_spectrum; % Element-wise operation, size N

% Generate road velocity profile (time domain)
road_velocity = ifft(final_spectrum, 'symmetric'); % Generate velocity signal
input_data = [time', road_velocity']; % Combine time and velocity as a 2-column matrix

% Plot road velocity profile
figure;
plot(time, road_velocity);
xlabel('Time (s)');
ylabel('Road Velocity (m/s)');
title('Generated Road Velocity Profile');

% Save road velocity profile for use in Simulink
save('road_velocity.mat', 'road_velocity', 'time');
