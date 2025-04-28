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


    
   
    % Creating time data
    STEP_SIMULATION = 0.0001; % Time step
    STEP_FREQUENCY = 0.1;   % Frequency resolution
    
    
    % Running model
    simOut = sim(model_name);
    
    

    % Loading Sprung mass acceleration data
    load("acceleration.mat","ans");
    time = [ans.Time(1):STEP_SIMULATION:ans.Time(end)]'; % Constant time step
    
    As = interp1(ans.Time, ans.Data, time); % Interpolated Sprung mass acceleration
    

    % Loading Tyre Force acceleration data
    load("tire_force.mat","ans");
    Ft = interp1(ans.Time, ans.Data, time); % Interpolated Tire Force



    As_filtered = lsim(tf([50 500], [1, 50, 1200]), As, time); %Acceleration sent through transfer function
    
    As_RMS = rms(As_filtered);
    Ft_RMS = rms(Ft);

    fprintf('Jr Value from Simulink Model (%s): %.3f m/s^2\n', model_name, As_RMS);
    fprintf('Jf Value from Simulink Model (%s): %.2f m/s^2\n', model_name, Ft_RMS);
    
    As_data = [time, As];
    
    % Process primary and secondary ride frequency ranges
    FFT_primary = FFT_JLR1(As_data, STEP_SIMULATION, STEP_FREQUENCY, time(1), time(end));
    FFT_secondary = FFT_JLR2(As_data, STEP_SIMULATION, STEP_FREQUENCY, time(1), time(end));
    FFT_msa_k(:,:)=FFT_primary;
    Lm1=length(FFT_primary);
    Lm2=length(FFT_secondary);
    FFT_msa_k(Lm1+1:Lm1+Lm2,1:2)=FFT_secondary;
    result = FFT_msa_k;

end


% clear all

% resultL0 = findFFT('mL0.slx');
% resultNL0 = findFFT('mNL0.slx');
% resultL1 = findFFT('mL1.slx');
% resultNL1 = findFFT('mNL1.slx');

resultNL1_hydraulic = findFFT('NL1_equib.slx');


% Plot the FFT result

% 
% figure(1);
% plot(resultL0(:,1),resultL0(:,2),'LineWidth',2, 'Color', 'red'); hold on;
% plot(resultL0_mechanical(:,1),resultL0_mechanical(:,2),'LineWidth',2, 'Color', 'magenta'); hold on;
% legend('L0 Hydraulic', 'L0 Mechanical');
% hold off;
% 
% 
% figure(2);
% plot(resultNL0(:,1),resultNL0(:,2),'LineWidth',2, 'Color', 'blue');hold on;
% plot(resultNL0_mechanical(:,1),resultNL0_mechanical(:,2),'LineWidth',2, 'Color', 'green'); hold on;
% legend('NL0 Hydraulic', 'NL0 Mechanical');
% hold off;

% figure(3);
% plot(resultL0(:,1),resultL0(:,2),'LineWidth',2, 'Color', 'blue'); hold on;
% plot(resultNL0(:,1),resultNL0(:,2),'LineWidth',2, 'Color','red' );hold on;
% plot(resultL1(:,1),resultL1(:,2),'LineWidth',2, 'Color', "#D95319"); hold on;
% plot(resultNL1(:,1),resultNL1(:,2),'LineWidth',2, 'Color', 'black' ); hold on;
% 
% xlabel('Frequency (Hz)')
% ylabel('J_r (ms^{-2})')
% 
% legend('L0', 'NL0', 'L1', 'NL1');
% % saveas(gcf,'Jr_performance_mechanical.png');
% hold off;


figure(6);
plot(resultNL1_hydraulic(:,1),resultNL1_hydraulic(:,2),'LineWidth',2, 'Color', 'blue'); hold on;
plot(resultNL1_hydraulic_parasitic(:,1),resultNL1_hydraulic(:,2),'LineWidth',2, 'Color','red' );hold on;

xlabel('Frequency (Hz)')
ylabel('J_r (ms^{-2})')

legend('L0', 'NL0', 'L1', 'NL1');
% saveas(gcf,'Jr_performance_mechanical.png');
hold off;
