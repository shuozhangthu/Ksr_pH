% Ion by ion model by Neilsen 213 with pH dependence from Wolthers 2012
%% global parameters
clear
% close all

load tang
load lorens
load tp

k_A10=8.2593;
v_A10=3;
a00=1.7;
epsilon0=-19.9875;
v_M10=3;
K_MB0=-7.6;
gamma0=-9.8268;
alpha0=1.41;

z0=[k_A10 v_A10  a00 epsilon0 K_MB0 gamma0 alpha0 v_M10]
lb=[1 0 a00 epsilon0 K_MB0 gamma0 alpha0 0];
ub=[9 5 a00 epsilon0 K_MB0 gamma0 alpha0 5];

 z = lsqcurvefit(@calculate_Rp_Ksr,z0,1,[log10(Rp_tang);ksr_tang;log10(R_lorens);Ksrca_lorens;log10(Rp_tp);ksr_tp]',lb,ub)

function output=calculate_Rp_Ksr(y,~)

K_AB=10^(-8.48);
K_MB=10^y(5);

a=6.4e-10;
b=a/2;
h=3.1e-10;
d=27100;
boltz=1.38065E-23;% Boltzmann constant
TC=25;
TK=273.15+TC;
gamma=10^y(6);
alpha=y(7);

a0=y(3);
epsilon=10^y(4);

%%
load tang
Omega=S_1_tang+1;

pH=8.3; %tang
M=5.1e-5;
A=0.0051;

theta=10^(8.6-pH);
phi=10^(10.33-pH);

k_A1=10^y(1);
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

v_A1=10^y(2);
v_A2=v_A1;
v_A=v_A1+v_A2;

v_M1=10^y(8);
v_M2=v_M1;
v_M=v_M1+v_M2;

for n=1:length(Omega)
    
    
    B=Omega(n)*K_AB/A;
    
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
        
        Kp8_high(n)=rx/M*A;
        
        u_BA=k_B*B*P_A-v_B*P_B*(1-x);
        u_BM=k_BM*B*P_M-v_BM*P_B*x;
        
        u_net=u_A+u_M+u_BA+u_BM;
        
        Omega_A=2+(v_A*exp(2*epsilon/boltz/TK)/k_B/B+v_A*exp(2*epsilon/boltz/TK)*v_B/k_A/A/k_B/B)/(1-v_A*v_B/k_A/A/k_B/B);
        Omega_B=2+(v_B*exp(2*epsilon/boltz/TK)/k_A/A+v_A*exp(2*epsilon/boltz/TK)*v_B/k_A/A/k_B/B)/(1-v_A*v_B/k_A/A/k_B/B);
        
        C=2/u_net*(k_A*A/Omega_A+k_B*B/Omega_B);
        
        rho_c=((C^2+4*C)^(0.5)-C)/2;
        
        v_st=rho_c*u_net*h/2;
        sigma=log(Omega(n));
        beta_st=v_st/sigma;
        Delta_G=1/8*h*a*b*gamma*gamma/h/boltz/TK/sigma;
        I=beta_st*h/h/a/b*(4*h*sigma/3.14/h/a/b)^0.5*exp(-Delta_G/boltz/TK);
        
        Rnet8_high(n) = 1.137*h*(I*v_st^2)^(1/3)*d; % Calcite growth rate
       
        y0=8*h*a*b*alpha/boltz/TK/sigma;
        Rnet8_low(n)=v_st*b*d/y0;

    end
    
    
    
end

%%
load lorens

Omega=Omega_lorens;

A=2.4e-3;
pH=7.4; %lorens
M=3.06e-7*0.24;%lorens
theta=10^(8.6-pH);
phi=10^(10.33-pH);

k_A1=10^y(1);
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


for n=1:length(Omega)
    
    
    B=Omega(n)*K_AB/A;
    
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
        sigma=log(Omega(n));        
        y0=8*h*a*b*alpha/boltz/TK/sigma;
        Rnet7_low(n)=v_st*b*d/y0;
                
    end
    
    
    
end




%%
load tp
Omega=Omega_tp;

pH=6.15;
A=3.8e-3;
M=1.46e-4;

theta=10^(8.6-pH);
phi=10^(10.33-pH);

k_A1=10^y(1);
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


for n=1:length(Omega)
    
    
    B=Omega(n)*K_AB/A;
    
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
        sigma=log(Omega(n));        
        y0=8*h*a*b*alpha/boltz/TK/sigma;
        Rnet6_low(n)=v_st*b*d/y0;
         
    end
    
    
    
end


 output=[log10(Rnet8_low) Kp8_high log10(Rnet7_low) Kp7_low log10(Rnet6_low) Kp6_low];

end
