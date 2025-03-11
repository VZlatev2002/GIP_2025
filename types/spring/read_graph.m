clear all


% load("NL1_a_100s.mat", 'Ms_acc');
load("types/spring/KL3_a.mat","Ms_acc")
time = Ms_acc.Time;
acc_data = Ms_acc.Data;

% load("NL1_H_100s.mat","H2631") H is the Jr
load("types/spring/KL3_H.mat","H2631")
H2631 = H2631.Data;
RMS= rms(H2631);
fprintf('Jr Value from Simulink Model: %.10f m/s^2\n', RMS);



%your_acceleration_data = [time, acc_data]; % Replace with your data

% Simulation Parameters
STEP_SIMULATION = 0.0001; % Time step
STEP_FREQUENCY = 0.1;   % Frequency resolution
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


function result = FFT_JLR1(a_sprung_ratley,STEP_SIMULATION,STEP_FREQUENCY,t1,t2)

[~,index_10s]=min(abs(a_sprung_ratley(:,1)-t1));
[~,index_40s]=min(abs(a_sprung_ratley(:,1)-t2));
x = a_sprung_ratley(index_10s:index_40s,2);

Fs=1/STEP_SIMULATION;

w_length1 = 20;%s 6 - flattopwin; 3 - hamming; 4 - blackman;
N_win1 = w_length1/STEP_SIMULATION;
window1 = flattopwin(N_win1);%Flat top window
noverlap1 = 0.85*N_win1;
f_range1 = [0:STEP_FREQUENCY:1.5];

[SxxWelch1,f1] = pwelch(x,window1,noverlap1,f_range1,Fs);
rmsDensity1 = SxxWelch1.^0.5;
result =[f1' rmsDensity1'];

end


function result = FFT_JLR2(a_sprung_ratley,STEP_SIMULATION,STEP_FREQUENCY,t1,t2)

[~,index_10s]=min(abs(a_sprung_ratley(:,1)-t1));
[~,index_40s]=min(abs(a_sprung_ratley(:,1)-t2));
x = a_sprung_ratley(index_10s:index_40s,2);

Fs=1/STEP_SIMULATION;

w_length2 = 6;%s 6 - flattopwin; 3 - hamming; 4 - blackman;
N_win2 = w_length2/STEP_SIMULATION;
window2 = flattopwin(N_win2);%Flat top window
noverlap2 = 0.85*N_win2;
f_range2 = [1.501:STEP_FREQUENCY:20];

[SxxWelch2,f2] = pwelch(x,window2,noverlap2,f_range2,Fs);
rmsDensity2 = SxxWelch2.^0.5;
result=[f2' rmsDensity2'];
end