% Load voltage and current data
voltageData = load('rlc-simple/voltage.mat');
currentData = load('rlc-simple/current.mat');

% Extract voltage and current values (assuming second row contains the data)
voltage = voltageData.voltage(2, :);
current = currentData.current(2, :);

figure;

% Plot the data
plot(-voltage, current, 'LineWidth', 1.5);

% Set axis labels and title
xlabel('Voltage (V)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Current (A)', 'FontSize', 12, 'FontWeight', 'bold');
title('Current-Voltage Characteristics', 'FontSize', 14, 'FontWeight', 'bold');

% Remove the box around the plot
ax = gca;
ax.Box = 'off'; % Turn off the box

% Optionally, you can set the line width of the axes
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;


% Set axis properties
set(gca, 'FontSize', 12);