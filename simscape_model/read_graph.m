clear all


load("NL5_a.mat", 'Ms_acc');

time = Ms_acc.Time;
acc_data = Ms_acc.Data;

load("NL5_H.mat","H2631")
H2631 = H2631.Data;
RMS= rms(H2631);
fprintf('Jr Value from Simulink Model: %.10f m/s^2\n', RMS);



%your_acceleration_data = [time, acc_data]; % Replace with your data

% Simulation Parameters
STEP_SIMULATION = 0.0001; % Time step
STEP_FREQUENCY = 0.01;   % Frequency resolution
t1 = 0;  % Start time for FFT
t2 = 100; % End time for FFT

tsimnew=[time(1):STEP_SIMULATION:time(end)]'; % to make a constant time step
msanew=interp1(time,acc_data,tsimnew); 

your_acceleration_data = [tsimnew,msanew];

% Process primary and secondary ride frequency ranges
FFT_primary = FFT_JLR1(your_acceleration_data, STEP_SIMULATION, STEP_FREQUENCY, t1, t2);
FFT_secondary = FFT_JLR2(your_acceleration_data, STEP_SIMULATION, STEP_FREQUENCY, t1, t2);


FFT_msa_k(:,:)=FFT_primary;

Lm1=length(FFT_primary);
Lm2=length(FFT_secondary);
FFT_msa_k(Lm1+1:Lm1+Lm2,1:2)=FFT_secondary;

figure;
% Plot the FFT result
%plot(FFT_result(:,1), FFT_result(:,2), 'LineWidth', 2);
plot(FFT_msa_k(:,1),FFT_msa_k(:,2),'LineWidth',2);hold on;
