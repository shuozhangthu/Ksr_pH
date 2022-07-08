 % wolthers ion by ion model supplemented by Neilsen 2012
clear
close all
set(0, 'DefaultAxesFontWeight', 'bold', ...
    'DefaultAxesFontSize', 16, ...
    'DefaultAxesFontAngle', 'normal', ... % Not sure the difference here
    'DefaultAxesFontWeight', 'bold', ... % Not sure the difference here
    'DefaultAxesTitleFontWeight', 'bold', ...
    'DefaultAxesTitleFontSizeMultiplier', 1) ;
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesLineWidth', 1.5)
set(0, 'DefaultLineMarkerSize', 6)

optimize_1D_jia;
%% global parameters
K_AB=10^(-8.48);
K_MB=10^z(5);

a=6.4e-10;
b=a/2;
h=3.1e-10;
d=27100;
boltz=1.38065E-23;% Boltzmann constant
TC=25;
TK=273.15+TC;
gamma=10^z(6);
alpha=z(7);

a0=z(3);
epsilon=10^z(4);

%% Tang
pH=8.3; %tang
M=5.1e-5;
A=0.0051;

theta=10^(8.6-pH);
phi=10^(10.33-pH);

super_sat=[0.01:0.1:2];
Omega8_low=10.^super_sat;

k_A1=10^z(1);
k_A2=k_A1;
k_B1=2*(1+theta)*k_A1/(1+phi);
k_B2=k_B1;
k_M1=0.3*k_A1;
k_M2=k_M1;
k_BM1=k_B1;
k_BM2=k_BM1;

k_A=k_A1+k_A2*theta;
k_B=k_B1+k_B2*phi;
k_M=k_M1+k_M2*theta;
k_BM=k_BM1+k_BM2*phi;


v_A1=10^z(2);
v_A2=v_A1;
v_A=v_A1+v_A2;

v_M1=10^z(8);
v_M2=v_M1;
v_M=v_M1+v_M2;

for n=1:length(Omega8_low)
    
    
    B=Omega8_low(n)*K_AB/A;
    
    P_A=0.5;
    P_M=0.01;
    P_B=(1-P_A-P_M)/(1+theta);
    x=0;
    
    
    
    for i=1:20
        
        v_B=exp(a0*x^2)*K_AB*k_A*k_B/v_A;
        v_BM=exp(a0*(1-x)^2)*K_MB*k_M*k_BM/v_M;
             
        [P_A,P_M,P_B]=solveP(k_A,k_B,k_M,k_BM,v_A,v_B,v_M,v_BM,P_A,P_M,A,B,M,theta);
        
        u_A=k_A*A*P_B-v_A*P_A;
        u_M=k_M*M*P_B-v_M*P_M;
        
        rx=u_M/u_A;
        x=rx/(1+rx);
        
        Kp8_low(n)=rx/M*A;
        
        u_BA=k_B*B*P_A-v_B*P_B*(1-x);
        u_BM=k_BM*B*P_M-v_BM*P_B*x;
        
        u_net=u_A+u_M+u_BA+u_BM;
        
        Omega_A=2+(v_A*exp(2*epsilon/boltz/TK)/k_B/B+v_A*exp(2*epsilon/boltz/TK)*v_B/k_A/A/k_B/B)/(1-v_A*v_B/k_A/A/k_B/B);
        Omega_B=2+(v_B*exp(2*epsilon/boltz/TK)/k_A/A+v_A*exp(2*epsilon/boltz/TK)*v_B/k_A/A/k_B/B)/(1-v_A*v_B/k_A/A/k_B/B);
        
        C=2/u_net*(k_A*A/Omega_A+k_B*B/Omega_B);
        
        rho_c=((C^2+4*C)^(0.5)-C)/2;
        v_st=rho_c*u_net*h/2;
        sigma=log(Omega8_low(n));
        y0=8*h*a*b*alpha/boltz/TK/sigma;
        Rnet8_low(n)=v_st*b*d/y0;
        
        RfRb=k_A*A*P_B/v_A/P_A;%Rf/Rb
        RpRb=RfRb-1;%Rp/Rb
        Rb_8(n)=Rnet8_low(n)/RpRb;
        Rf_8(n)=Rb_8(n)*RfRb;
    end
    PB1_8(n)=P_B;
    v_st_8(n)=v_st;
    k_B_8(n)=k_B;
    PA_8(n)=P_A;
    PM_8(n)=P_M;
end

%%
%Lorens
A=2.4e-3;
pH=7.4; %lorens
M=3.06e-7*0.24;%lorens

theta=10^(8.6-pH);
phi=10^(10.33-pH);

super_sat=[0.01:0.1:2];
Omega7_low=10.^super_sat;

k_A1=10^z(1);
k_A2=k_A1;
k_B1=2*(1+theta)*k_A1/(1+phi);
k_B2=k_B1;
k_M1=0.3*k_A1;
k_M2=k_M1;
k_BM1=k_B1;
k_BM2=k_BM1;

k_A=k_A1+k_A2*theta;
k_B=k_B1+k_B2*phi;
k_M=k_M1+k_M2*theta;
k_BM=k_BM1+k_BM2*phi;


for n=1:length(Omega7_low)
    
    
    B=Omega7_low(n)*K_AB/A;
    
    P_A=0.5;
    P_M=0.01;
    P_B=(1-P_A-P_M)/(1+theta);
    x=0;
    
    
    
    for i=1:20
        
        v_B=exp(a0*x^2)*K_AB*k_A*k_B/v_A;
        v_BM=exp(a0*(1-x)^2)*K_MB*k_M*k_BM/v_M;
      
        [P_A,P_M,P_B]=solveP(k_A,k_B,k_M,k_BM,v_A,v_B,v_M,v_BM,P_A,P_M,A,B,M,theta);
        
        u_A=k_A*A*P_B-v_A*P_A;
        u_M=k_M*M*P_B-v_M*P_M;
        
        rx=u_M/u_A;
        x=rx/(1+rx);
        
        Kp7_low(n)=rx/M*A;
        
        u_BA=k_B*B*P_A-v_B*P_B*(1-x);
        u_BM=k_BM*B*P_M-v_BM*P_B*x;
        
        u_net=u_A+u_M+u_BA+u_BM;
        
        Omega_A=2+(v_A*exp(2*epsilon/boltz/TK)/k_B/B+v_A*exp(2*epsilon/boltz/TK)*v_B/k_A/A/k_B/B)/(1-v_A*v_B/k_A/A/k_B/B);
        Omega_B=2+(v_B*exp(2*epsilon/boltz/TK)/k_A/A+v_A*exp(2*epsilon/boltz/TK)*v_B/k_A/A/k_B/B)/(1-v_A*v_B/k_A/A/k_B/B);
        
        C=2/u_net*(k_A*A/Omega_A+k_B*B/Omega_B);
        
        rho_c=((C^2+4*C)^(0.5)-C)/2;
        v_st=rho_c*u_net*h/2;
        sigma=log(Omega7_low(n));
        y0=8*h*a*b*alpha/boltz/TK/sigma;
        Rnet7_low(n)=v_st*b*d/y0;
        RfRb=k_A*A*P_B/v_A/P_A;%Rf/Rb
        RpRb=RfRb-1;%Rp/Rb
        Rb_7(n)=Rnet7_low(n)/RpRb;
        Rf_7(n)=Rb_7(n)*RfRb;
    end
    
    PB1_7(n)=P_B;
    v_st_7(n)=v_st;
    k_B_7(n)=k_B;
    PA_7(n)=P_A;
    PM_7(n)=P_M;
end



%%
%T&P
pH=6.15;
A=3.8e-3;
M=1.46e-4;

theta=10^(8.6-pH);
phi=10^(10.33-pH);
super_sat=[0.01:0.1:2];
Omega6_low=10.^super_sat;
k_A1=10^z(1);
k_A2=k_A1;
k_B1=2*(1+theta)*k_A1/(1+phi);
k_B2=k_B1;
k_M1=0.3*k_A1;
k_M2=k_M1;
k_BM1=k_B1;
k_BM2=k_BM1;

k_A=k_A1+k_A2*theta;
k_B=k_B1+k_B2*phi;
k_M=k_M1+k_M2*theta;
k_BM=k_BM1+k_BM2*phi;



for n=1:length(Omega6_low)
    
    
    B=Omega6_low(n)*K_AB/A;
    
    P_A=0.5;
    P_M=0.01;
    P_B=(1-P_A-P_M)/(1+theta);
    x=0;
    
    
    
    for i=1:20
        v_B=exp(a0*x^2)*K_AB*k_A*k_B/v_A;
        v_BM=exp(a0*(1-x)^2)*K_MB*k_M*k_BM/v_M;  
        [P_A,P_M,P_B]=solveP(k_A,k_B,k_M,k_BM,v_A,v_B,v_M,v_BM,P_A,P_M,A,B,M,theta);
        
        u_A=k_A*A*P_B-v_A*P_A;
        u_M=k_M*M*P_B-v_M*P_M;
        
        rx=u_M/u_A;
        x=rx/(1+rx);
        
        Kp6_low(n)=rx/M*A;
        
        u_BA=k_B*B*P_A-v_B*P_B*(1-x);
        u_BM=k_BM*B*P_M-v_BM*P_B*x;
        
        u_net=u_A+u_M+u_BA+u_BM;
        
        Omega_A=2+(v_A*exp(2*epsilon/boltz/TK)/k_B/B+v_A*exp(2*epsilon/boltz/TK)*v_B/k_A/A/k_B/B)/(1-v_A*v_B/k_A/A/k_B/B);
        Omega_B=2+(v_B*exp(2*epsilon/boltz/TK)/k_A/A+v_A*exp(2*epsilon/boltz/TK)*v_B/k_A/A/k_B/B)/(1-v_A*v_B/k_A/A/k_B/B);
        
        C=2/u_net*(k_A*A/Omega_A+k_B*B/Omega_B);
        
        rho_c=((C^2+4*C)^(0.5)-C)/2;
        v_st=rho_c*u_net*h/2;
        sigma=log(Omega6_low(n));
        y0=8*h*a*b*alpha/boltz/TK/sigma;
        Rnet6_low(n)=v_st*b*d/y0;
        
                   
        RfRb=k_A*A*P_B/v_A/P_A;%Rf/Rb
        RpRb=RfRb-1;%Rp/Rb
        Rb_6(n)=Rnet6_low(n)/RpRb;
        Rf_6(n)=Rb_6(n)*RfRb;

    end
    PB1_6(n)=P_B;
    PA_6(n)=P_A;
    PM_6(n)=P_M;

    v_st_6(n)=v_st;
    k_B_6(n)=k_B;
end
%% Experimental data
load tang
load lorens
load tp
%% plot figures

% 
figure
set(gcf,'Position',[100, 0, 700,600])
hold on
scatter(S_1_tang,Rp_tang,60,[0.9216,0.3020,0.2941],'s','MarkerFaceColor',[1,0.4745,0.4745],'linewidth',1.5)
scatter(Omega_lorens-1,R_lorens,60,[0.9412,0.5765,0.1686],'o','MarkerFaceColor',[1,0.7451,0.4627],'linewidth',1.5)
scatter(Omega_tp-1,Rp_tp,60,[0.4157,0.6902,0.2980],'^','MarkerFaceColor','(0.7294,0.8627,0.3451)','linewidth',1.5)
plot(Omega8_low-1, Rnet8_low,'--','color',[0.9216,0.3020,0.2941],'linewidth',2)
plot(Omega7_low-1, Rnet7_low,'-.','color',[0.9412,0.5765,0.1686],'linewidth',2)
plot(Omega6_low-1, Rnet6_low,'color',[0.4157,0.6902,0.2980],'linewidth',2)
axis([0.1 100 1*10^(-12) 1*10^(-3)])
xlabel('\Omega-1','fontsize',18,'interpreter','tex')
ylabel('R_p (mol/m^2/s)','fontsize',18,'interpreter','tex')
set(gca,'xscale','log')
set(gca,'yscale','log')
box on
ax = gca;
ax.LineWidth = 1.5;
legend('Tang et al. (2008), pH= 8.3','Lorens (1981), pH= 7.4','T&P (1996), pH= 6.15','model, pH= 8.3, spiral','model, pH= 7.4, spiral','model, pH=6.15, spiral','location','best')
ylim([1e-10,1e-4]);
xlim([0.01,100]);
print('Rp_S.jpeg','-djpeg','-r1200');
set(findobj(gcf,'type','axes'),'FontName','Times New Roman')



figure
set(gcf,'Position',[100, 0, 700,600])
hold on
scatter(S_1_tang,ksr_tang,60,[0.9216,0.3020,0.2941],'s','MarkerFaceColor',[1,0.4745,0.4745],'linewidth',1.5)
scatter(Omega_lorens-1,Ksrca_lorens,60,[0.9412,0.5765,0.1686],'o','MarkerFaceColor',[1,0.7451,0.4627],'linewidth',1.5)
scatter(Omega_tp-1,ksr_tp,60,[0.4157,0.6902,0.2980],'^','MarkerFaceColor','(0.7294,0.8627,0.3451)','linewidth',1.5)
plot(Omega8_low-1,Kp8_low,'--','color',[0.9216,0.3020,0.2941],'linewidth',2)
plot(Omega7_low-1, Kp7_low,'-.', 'color',[0.9412,0.5765,0.1686],'linewidth',2)
plot(Omega6_low-1, Kp6_low, 'color',[0.4157,0.6902,0.2980],'linewidth',2)
xlabel('\Omega-1','fontsize',18,'interpreter','tex')
ylabel('K_{Sr}','fontsize',18,'interpreter','tex')
set(gca,'xscale','log')
box on
ax = gca;
ax.LineWidth = 1.5;
legend('Tang et al. (2008), pH= 8.3','Lorens (1981), pH= 7.4','T&P (1996), pH= 6.15','model, pH= 8.3, spiral','model, pH= 7.4, spiral','model, pH=6.15, spiral','location','best')
xlim([1e-2,1e2]);
ylim([0,0.25]);
set(findobj(gcf,'type','axes'),'FontName','Times New Roman')

print('Ksr_S.jpeg','-djpeg','-r1200');


figure
set(gcf,'Position',[100, 0, 700,600])
hold on
scatter(Rp_tang,ksr_tang,60,[0.9216,0.3020,0.2941],'s','MarkerFaceColor',[1,0.4745,0.4745],'linewidth',1.5)
scatter(R_lorens,Ksrca_lorens,60,[0.9412,0.5765,0.1686],'o','MarkerFaceColor',[1,0.7451,0.4627],'linewidth',1.5)
scatter(Rp_tp,ksr_tp,60,[0.4157,0.6902,0.2980],'^','MarkerFaceColor','(0.7294,0.8627,0.3451)','linewidth',1.5)
plot(Rnet8_low, Kp8_low,'--','color',[0.9216,0.3020,0.2941],'linewidth',2)
plot(Rnet7_low, Kp7_low,'-.', 'color',[0.9412,0.5765,0.1686],'linewidth',2)
plot(Rnet6_low, Kp6_low,'color',[0.4157,0.6902,0.2980],'linewidth',2)
xlabel('R_p (mol/m^2/s)','fontsize',18,'interpreter','tex')
ylabel('K_{Sr}','fontsize',18,'interpreter','tex')
set(gca,'xscale','log')
box on
ax = gca;
ax.LineWidth = 1.5;
legend('Tang et al. (2008), pH= 8.3','Lorens (1981), pH= 7.4','T&P (1996), pH= 6.15','model, pH= 8.3, spiral','model, pH= 7.4, spiral','model, pH=6.15, spiral','location','best')
xlim([1e-10,1e-4]);
ylim([0,0.25]);
set(findobj(gcf,'type','axes'),'FontName','Times New Roman')

print('Ksr_Rp.jpeg','-djpeg','-r1200');

% 
% 
figure
set(gcf,'Position',[100, 0, 800,1000])
subplot(2,2,1)
hold on
PB2_8=PB1_8.*10^(8.6-8.3);
P_B8=PB2_8+PB1_8;
plot(Omega8_low-1,PB1_8,'--','color',[0.9216,0.3020,0.2941],'linewidth',1.5)
PB2_7=PB1_7.*10^(8.6-7.4);
P_B7=PB2_7+PB1_7;
plot(Omega7_low-1,PB1_7,'-.','color',[0.9412,0.5765,0.1686],'linewidth',1.5)
PB2_6=PB1_6.*10^(8.6-6.15);
P_B6=PB2_6+PB1_6;
plot(Omega6_low-1,PB1_6,'color',[0.4157,0.6902,0.2980],'linewidth',1.5)
box on
ax = gca;
ax.LineWidth = 1;
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('\Omega-1','interpreter','tex')
ylabel('P_{B1}','interpreter','tex')
set(gca,'FontSize',10) 
xlabel('\Omega-1','fontsize',12,'interpreter','tex')
ylabel('P_{B1}','fontsize',12,'interpreter','tex')
ylim([0.0001,0.1]);
xlim([0.01,100]);
title(['(a)']);
set(findobj(gcf,'type','axes'),'FontName','Times New Roman')

subplot(2,2,2)
plot(Omega8_low-1,PA_8,'--','color',[0.9216,0.3020,0.2941],'linewidth',1.5)
hold on
plot(Omega7_low-1,PA_7,'-.', 'color',[0.9412,0.5765,0.1686],'linewidth',1.5)
plot(Omega6_low-1,PA_6,'color',[0.4157,0.6902,0.2980],'linewidth',1.5)
set(gca,'xscale','log')
ylim([0.1,10]);
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'FontSize',10) 
xlabel('\Omega-1','fontsize',12,'interpreter','tex')
ylabel('P_{A}','fontsize',12,'interpreter','tex')
box on
ax = gca;
ax.LineWidth = 1;
title(['(b)'],'Fontsize',10);
set(findobj(gcf,'type','axes'),'FontName','Times New Roman')

subplot(2,2,3)
hold on
plot(Omega8_low-1,PM_8,'--','color',[0.9216,0.3020,0.2941],'linewidth',1.5)
plot(Omega7_low-1,PM_7,'-.', 'color',[0.9412,0.5765,0.1686],'linewidth',1.5)
plot(Omega6_low-1,PM_6,'color',[0.4157,0.6902,0.2980],'linewidth',1.5)
set(gca,'xscale','log')
ylim([0.000001,0.1]);
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'FontSize',10) 
xlabel('\Omega-1','fontsize',12,'interpreter','tex')
ylabel(' P_{M}','fontsize',12,'interpreter','tex')
box on
ax = gca;
ax.LineWidth = 1;
title(['(c)'],'Fontsize',10);
title(['(c)']);
set(findobj(gcf,'type','axes'),'FontName','Times New Roman')


subplot(2,2,4)
plot(Omega8_low-1,Rb_8./Rf_8,'--','color',[0.9216,0.3020,0.2941],'linewidth',1.5)
hold on
plot(Omega7_low-1,Rb_7./Rf_7,'-.', 'color',[0.9412,0.5765,0.1686],'linewidth',1.5)
plot(Omega6_low-1,Rb_6./Rf_6,'color',[0.4157,0.6902,0.2980],'linewidth',1.5)
set(gca,'xscale','log')
set(gca,'FontSize',10) 
xlabel('\Omega-1','fontsize',12,'interpreter','tex')
ylabel('R_{b Ca^{2+}}/R_{f Ca^{2+}}','fontsize',12,'interpreter','tex')
box on
box on
ax = gca;
ax.LineWidth = 1;
title(['(d)'],'Fontsize',10);
legend('Tang et al. (2008), pH = 8.3','Lorens (1981), pH = 7.4','T&P (1996), pH = 6.15','location','best')
set(findobj(gcf,'type','axes'),'FontName','Times New Roman')

print('P.jpeg','-djpeg','-r1200');


