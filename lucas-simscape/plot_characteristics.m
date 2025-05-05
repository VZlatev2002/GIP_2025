VDD = sim("VDD.slx");



% VDD_blowoff_valve = sim("VDD_blowoff_valve.slx");
% ZPF = sim("ZPF.slx");
% Passive_damper = sim("Passive_damper.slx");


% % COMPONENT TEST PLOT
% figure(3);
% component_test = sim("single_component_testing.slx");
% plot(component_test.VDD_velocity.Data, component_test.VDD_force.Data,'LineWidth',2, 'Color', 'blue'); hold on;
% xlabel('Velocity'); hold on
% grid on
% ylabel('Force') 


figure(1);
plot(VDD.VDD_displacement.Data, VDD.VDD_force.Data,'LineWidth',2, 'Color', '#008000'); hold on;
xlabel('Displacement (m)') 
ylabel('Force (N)') 

figure(4);
plot(VDD.VDD_velocity.Data, VDD.VDD_force.Data,'LineWidth',2, 'Color', '#008000'); hold on;
xlabel('Velocity (m/s)') 
ylabel('Force (N)') 

