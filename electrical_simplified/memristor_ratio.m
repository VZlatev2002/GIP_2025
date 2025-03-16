% Define R_A (high resistance) and R_B (low resistance)
R_A = 17; % Ω (example)
R_B = 7;   % Ω (example)
simOut = sim('new_electrical.slx');

% Extract xi from simulation output
xiData = simOut.logsout.get('xi').Values;
time = xiData.Time;
xi = xiData.Data;

% Compute resistance R(t)
R = xi .* R_A + (1 - xi) .* R_B;

% Plot resistance over time
plot(time, xi);
xlabel('Time (s)');
ylabel('Memristor Resistance (Ω)');
title('Computed Memristor Resistance Over Time');
grid on;
