% Load voltage and current data
voltageData = load('electrical_simplified/voltage.mat');
currentData = load('electrical_simplified/current.mat');

% Extract voltage and current values (assuming second row contains the data)
voltage = voltageData.voltage(2, :);
current = currentData.current(2, :);

% Create a new figure
figure;

% --- Subplot 1: Scatter Plot ---
subplot(1, 2, 1);  % 1 row, 2 columns, use the first cell
scatter(-voltage, current, 10, 'b', 'filled');  % Negative voltage if desired
xlabel('Voltage (V)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Current (A)', 'FontSize', 12, 'FontWeight', 'bold');
title('Scatter Plot', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);

% --- Subplot 2: Continuous Line Plot ---
subplot(1, 2, 2);  % Still the same figure, second cell
plot(-voltage, current, 'LineWidth', 1.5);
xlabel('Voltage (V)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Current (A)', 'FontSize', 12, 'FontWeight', 'bold');
title('Line Plot', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);

% Save the entire figure as an image
saveas(gcf, 'IV_Curve_TwoSubplots.png');

disp('Scatter and continuous line plots created on two separate axes, saved as IV_Curve_TwoSubplots.png');
