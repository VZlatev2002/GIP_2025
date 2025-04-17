%% Clean workspace
clear all
clc

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

function result = findFFT(model_name)

    simOut = sim(model_name);
    load("rlc_simple_two/a.mat","ans");
    time = ans.Time;
    acc_data = ans.Data;
    
    
    % load("NL1_H_100s.mat","H2631") H is the Jr
    load("rlc_simple_two/h.mat","H2631")
    H2631 = H2631.Data;
    RMS= rms(H2631);
    fprintf('Jr Value from Simulink Model: %.10f m/s^2\n', RMS);

    % Simulation Parameters
    STEP_SIMULATION = 0.0001; % Time step
    STEP_FREQUENCY = 0.1;   % Frequency resolution
    t1 = 0;  % Start time for FFT
    t2 = 100; % End time for FFT
    
    tsimnew=[time(1):STEP_SIMULATION:time(end)]'; % to make a constant time step
    msanew=interp1(time,acc_data,tsimnew); 
    
    filtered_acceleration = lsim(tf([50 500], [1, 50, 1200]), msanew, tsimnew);
    
    RMS = rms(filtered_acceleration); 
 %   fprintf('Jr Value from Simulink Model (%s): %.2f m/s^2\n', model_name, RMS);
    your_acceleration_data = [tsimnew,msanew];
    
    % Process primary and secondary ride frequency ranges
    FFT_primary = FFT_JLR1(your_acceleration_data, STEP_SIMULATION, STEP_FREQUENCY, t1, t2);
    FFT_secondary = FFT_JLR2(your_acceleration_data, STEP_SIMULATION, STEP_FREQUENCY, t1, t2);
    
    
    FFT_msa_k(:,:)=FFT_primary;
    
    Lm1=length(FFT_primary);
    Lm2=length(FFT_secondary);
    FFT_msa_k(Lm1+1:Lm1+Lm2,1:2)=FFT_secondary;
    result = FFT_msa_k;

end


clear all`


resultNL0 = findFFT('rlc_full_two.slx');


figure(2);
plot(resultNL0(:,1),resultNL0(:,2),'LineWidth',2, 'Color', 'blue');hold on;
legend('NL1 Electrical', 'NL1 Electrical')
hold off;