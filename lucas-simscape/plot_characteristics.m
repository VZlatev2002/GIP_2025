VDD = sim("VDD.slx");
ZPF = sim("ZPF.slx");
Passive_damper = sim("Passive_damper.slx");
figure(1);



% figure(2);
plot(ZPF.ZPF_velocity.Data(300:1000), ZPF.ZPF_force.Data(300:1000),'LineWidth',2, 'Color', 'red'); hold on;
%set(gca, 'XTick', [], 'YTick', []);
xlabel('Velocity') 
ylabel('Force') 

% scatter(Passive_damper.Passive_velocity.Data, Passive_damper.Passive_force.Data,'LineWidth',2, 'Color', 'blue'); hold on;
plot(VDD.VDD_velocity.Data(300:30000), VDD.VDD_force.Data(300:30000),'LineWidth',2, 'Color', 'blue'); hold on;
xlabel('Velocity'); hold on
ylabel('Force') 
%set(gca, 'XTick', [], 'YTick', []);
% saveas(gcf,'Hydraulic_VDD_velocity-force.png')
legend('ZPF', 'VDD', 'Location','southeast')
%saveas(gcf,'ZPF_VDD_velocity-force.png')


figure(3);

% figure(4);
plot(ZPF.ZPF_displacement.Data(300:1000), ZPF.ZPF_force.Data(300:1000),'LineWidth',2, 'Color', 'red'); hold on;
%set(gca, 'XTick', [], 'YTick', []);
xlabel('Displacement') 
ylabel('Force') 

% scatter(Passive_damper.Passive_displacement.Data, Passive_damper.Passive_force.Data,'LineWidth',2, 'Color', 'blue'); hold on;
plot(VDD.VDD_displacement.Data(300:30000), VDD.VDD_force.Data(300:30000),'LineWidth',2, 'Color', 'blue'); hold on;
%set(gca, 'XTick', [], 'YTick', []);
xlabel('Displacement') 
ylabel('Force') 
% saveas(gcf,'Hydraulic_VDD_displacement-force.png')
legend('ZPF', 'VDD')
%saveas(gcf,'ZPF_VDD_displacement-force.png')
