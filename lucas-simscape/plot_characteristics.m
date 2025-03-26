VDD = sim("VDD.slx");
ZPF = sim("ZPF.slx");
Passive_damper = sim("Passive_damper.slx");
figure(1);
scatter(Passive_damper.Passive_velocity.Data, Passive_damper.Passive_force.Data,'LineWidth',2, 'Color', 'blue'); hold on;
scatter(VDD.VDD_velocity.Data, VDD.VDD_force.Data,'LineWidth',2, 'Color', 'blue'); hold on;
scatter(ZPF.ZPF_velocity.Data, ZPF.ZPF_force.Data,'LineWidth',2, 'Color', 'blue'); hold off;

figure(3);
scatter(Passive_damper.Passive_displacement.Data, Passive_damper.Passive_force.Data,'LineWidth',2, 'Color', 'blue'); hold on;
scatter(VDD.VDD_displacement.Data, VDD.VDD_force.Data,'LineWidth',2, 'Color', 'blue'); hold on;
scatter(ZPF.ZPF_displacement.Data, ZPF.ZPF_force.Data,'LineWidth',2, 'Color', 'blue'); hold off;