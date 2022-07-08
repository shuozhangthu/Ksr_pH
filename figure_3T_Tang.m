clear
run_5;
load tang_all
figure
set(gcf,'Position',[100, 0, 700,600])
hold on
scatter(omega_tang5-1, ksr_tang5,60,[0.9412,0.5765,0.1686],'o','MarkerFaceColor',[1,0.7451,0.4627],'linewidth',1.5)
scatter(omega_tang25-1, ksr_tang25,60,[0.9216,0.3020,0.2941],'s','MarkerFaceColor',[1,0.4745,0.4745],'linewidth',1.5)
scatter(omega_tang40-1, ksr_tang40,60,[0.4157,0.6902,0.2980],'^','MarkerFaceColor','(0.7294,0.8627,0.3451)','linewidth',1.5)
set(gca,'ColorOrderIndex',1)

plot(Omega8_low-1, Kp8_low,'-.','color',[0.9412,0.5765,0.1686],'linewidth',2)
set(gca,'xscale','log')
box on
ax = gca;
ax.LineWidth = 1.5;

clear
run_25;
hold on
plot(Omega8_low-1, Kp8_low,'--','color',[0.9216,0.3020,0.2941],'linewidth',2)
set(gca,'xscale','log')
box on
ax = gca;
ax.LineWidth = 1.5;


clear
run_40;
load tang_all
N = 3;
hold on 
plot(Omega8_low-1, Kp8_low,'color',[0.4157,0.6902,0.2980],'linewidth',2)
xlabel('\Omega-1','fontsize',18,'interpreter','tex')
ylabel('K_{Sr}','fontsize',18,'interpreter','tex')
set(gca,'xscale','log')
box on
ax = gca;
ax.LineWidth = 1.5;
set(findobj(gcf,'type','axes'),'FontName','Times New Roman')
legend('T = 5°C','T = 25°C','T = 40°C','2D nucl., T = 5°C','2D nucl., T = 25°C','2D nucl., T = 40°C','location','best')
print('Ksr_S_temp.jpeg','-djpeg','-r1200');

%%
figure
run_5;
set(gcf,'Position',[100, 0, 700,600])
hold on
scatter(omega_tang5-1,Rp_tang5,60,[0.9412,0.5765,0.1686],'o','MarkerFaceColor',[1,0.7451,0.4627],'linewidth',1.5)
scatter(omega_tang25-1,Rp_tang25,60,[0.9216,0.3020,0.2941],'s','MarkerFaceColor',[1,0.4745,0.4745],'linewidth',1.5)
scatter(omega_tang40-1,Rp_tang40,60,[0.4157,0.6902,0.2980],'^','MarkerFaceColor','(0.7294,0.8627,0.3451)','linewidth',1.5)
plot(Omega8_low-1,Rnet8_low,'-.','color',[0.9412,0.5765,0.1686],'linewidth',2)
xlabel('R_p (mol/m^2/s)')
ylabel('K_{Sr}')
set(gca,'xscale','log')
set(gca,'yscale','log')
box on
ax = gca;
ax.LineWidth = 1.5;

clear
run_25;
hold on
plot(Omega8_low-1,Rnet8_low,'--','color',[0.9216,0.3020,0.2941],'linewidth',2)
xlabel('R_p (mol/m^2/s)')
ylabel('K_{Sr}')
set(gca,'xscale','log')
box on
ax = gca;
ax.LineWidth = 1.5;


clear
run_40;
hold on 
plot(Omega8_low-1,Rnet8_low,'color',[0.4157,0.6902,0.2980],'linewidth',2)
xlabel('\Omega-1','fontsize',18,'interpreter','tex')
ylabel('R_p (mol/m^2/s)','fontsize',18,'interpreter','tex')
set(gca,'xscale','log')
box on
ax = gca;
ax.LineWidth = 1.5;
legend('T = 5°C','T = 25°C','T = 40°C','2D nucl., T = 5°C','2D nucl., T = 25°C','2D nucl., T = 40°C','location','best')
ylim([1e-10,1e-4]);
set(findobj(gcf,'type','axes'),'FontName','Times New Roman')
print('Rp_S_temp.jpeg','-djpeg','-r1200');

%%
figure
set(gcf,'Position',[100, 0, 700,600])
clear
run_5;
hold on
scatter(Rp_tang5,ksr_tang5,60,[0.9412,0.5765,0.1686],'o','MarkerFaceColor',[1,0.7451,0.4627],'linewidth',1.5)
scatter(Rp_tang25,ksr_tang25,60,[0.9216,0.3020,0.2941],'s','MarkerFaceColor',[1,0.4745,0.4745],'linewidth',1.5)
scatter(Rp_tang40,ksr_tang40,60,[0.4157,0.6902,0.2980],'^','MarkerFaceColor','(0.7294,0.8627,0.3451)','linewidth',1.5)
plot(Rnet8_low, Kp8_low,'-.','color',[0.9412,0.5765,0.1686],'linewidth',2)
xlabel('R_p (mol/m^2/s)')
ylabel('K_{Sr}')
set(gca,'xscale','log')
box on
ax = gca;
ax.LineWidth = 1.5;

clear
run_25;
hold on
plot(Rnet8_low, Kp8_low,'--','color',[0.9216,0.3020,0.2941],'linewidth',2)
xlabel('R_p (mol/m^2/s)')
ylabel('K_{Sr}')
set(gca,'xscale','log')
box on
ax = gca;
ax.LineWidth = 1.5;

clear
run_40;
hold on 
plot(Rnet8_low, Kp8_low,'color',[0.4157,0.6902,0.2980],'linewidth',2)
xlabel('R_p (mol/m^2/s)','fontsize',18,'interpreter','tex')
ylabel('K_{Sr}','fontsize',18,'interpreter','tex')
set(gca,'xscale','log')
box on
ax = gca;
ax.LineWidth = 1.5;
legend('T = 5°C','T = 25°C','T = 40°C','2D nucl., T = 5°C','2D nucl., T = 25°C','2D nucl., T = 40°C','location','best')
xlim([1e-9,1e-5]);
set(findobj(gcf,'type','axes'),'FontName','Times New Roman')
print('Ksr_Rp_temp.jpeg','-djpeg','-r1200');


