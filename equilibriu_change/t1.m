% Assuming 'out' is the output from the sim function
% Access the logged data
piston_disp_series = simlog.Double_Acting_Actuator_IL1.p_out.series;

% Extract time and data
time = piston_disp_series.time;
values = piston_disp_series.values;

% Plot the data
figure;
plot(time, values);
title('Piston Displacement');
xlabel('Time [s]');
ylabel('Displacement [m]');
grid on;
length(values)
mean(values)
std(values)
