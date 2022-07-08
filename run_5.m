% wolthers ion by ion model supplemented by Neilsen 2012
% clear
% close all
clear
set(0, 'DefaultAxesFontWeight', 'bold', ...
    'DefaultAxesFontSize', 14, ...
    'DefaultAxesFontAngle', 'normal', ... % Not sure the difference here
    'DefaultAxesFontWeight', 'bold', ... % Not sure the difference here
    'DefaultAxesTitleFontWeight', 'bold', ...
    'DefaultAxesTitleFontSizeMultiplier', 1) ;
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesLineWidth', 1.5)
set(0, 'DefaultLineMarkerSize', 6)

 optimize_5;
  z=[   10.0000   4.7367    4.9587];
 %% global parameters
K_AB=10^(-8.39);
K_MB=2.9627E-08;
a=6.4e-10;
b=a/2;
h=3.1e-10;
d=27100;
boltz=1.38065E-23;% Boltzmann constant
TC=5;
TK=273.15+TC;
alpha=1.41;
a0=1.7;
epsilon=10^(-19.9875);
gamma=1.49e-10;

%%
pH=8.99; %tang
M=5.275e-5;
A=0.005229;

theta=10^(8.6-pH);
phi=10^(10.55-pH);

super_sat=[0.3010:0.1:1.5];
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
v_M1=10^z(3);
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
        
        v_st=rho_c*u_net*b;
        sigma=log(Omega8_low(n));       
        y0=8*h*a*b*alpha/boltz/TK/sigma;   
          beta_st=v_st/sigma;
        Delta_G=1/8*h*a*b*gamma*gamma/h/boltz/TK/sigma;
        I=beta_st*h/h/a/b*(4*h*sigma/3.14/h/a/b)^0.5*exp(-Delta_G/boltz/TK);
        
        Rnet8_low(n) = 1.137*h*(I*v_st^2)^(1/3); % Calcite growth rate
     
    end
    
    
    
end

%% Experimental data
 load tang_all
