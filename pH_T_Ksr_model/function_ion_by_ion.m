%function [Rp,Ksr]=function_ion_by_ion(pH,omega,T,F)
%F=1,spiral ; F=2,2D
% the activities of Ca and Sr are consistent with that in seawater.
function [Rp,Ksr]=function_ion_by_ion(pH,omega,T,F)%% global parameters
T=25;
pH=8;
Omega=10;
TK=273.15+T;

K_AB=10^(-171.9065-0.077993*TK+2839.3/TK+71.595*log10(TK));

ksr_TK=0.1329*exp(5000/8.315*(1/298.15-1/TK)); 
lnKsrco3_TK=log(K_AB)-log(ksr_TK);
K_MB=exp(lnKsrco3_TK);

a=6.4e-10;
b=a/2;
h=3.1e-10;
d=27100;
boltz=1.38065E-23;% Boltzmann constant

gamma=10^(-9.7570);
alpha=1.41;

a0=1.7;
epsilon=10^(-19.9875);

A=0.002;
M=A*8.54/1000;
theta=10^(8.6-pH);

logKhco3_TK=-107.8871-0.03252894*TK+5151.79/TK+38.92561*log10(TK)-563713.9/TK/TK;
phi=10^(-logKhco3_TK-pH);

 if F==1
k_A1=10^(5.9653);%spiral
 else 
k_A1=884000000000*exp(-1280/TK);%2D
 end
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

 if F==1
v_A1=10^1.9681;%spiral
 else
v_A1=1598000*exp(-930.7/TK);%2D
 end
v_A2=v_A1;
v_A=v_A1+v_A2;

 if F==1
v_M1=10^2.0883;%spiral
 else
v_M1=16030000000000*exp(-5396/TK);%2D
 end

v_M2=v_M1;
v_M=v_M1+v_M2;



B=omega*K_AB/A;

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
    
    Ksr=rx/M*A;
    
    u_BA=k_B*B*P_A-v_B*P_B*(1-x);
    u_BM=k_BM*B*P_M-v_BM*P_B*x;
    
    u_net=u_A+u_M+u_BA+u_BM;
    
    Omega_A=2+(v_A*exp(2*epsilon/boltz/TK)/k_B/B+v_A*exp(2*epsilon/boltz/TK)*v_B/k_A/A/k_B/B)/(1-v_A*v_B/k_A/A/k_B/B);
    Omega_B=2+(v_B*exp(2*epsilon/boltz/TK)/k_A/A+v_A*exp(2*epsilon/boltz/TK)*v_B/k_A/A/k_B/B)/(1-v_A*v_B/k_A/A/k_B/B);
    
    C=2/u_net*(k_A*A/Omega_A+k_B*B/Omega_B);
    
    rho_c=((C^2+4*C)^(0.5)-C)/2;
    
    v_st=rho_c*u_net*h/2;
    sigma=log(omega);

if F==1
%spiral
    y0=8*h*a*b*alpha/boltz/TK/sigma;
    Rp=v_st*b*d/y0;
else
    %2D
    beta_st=v_st/sigma;
    Delta_G=1/8*h*a*b*gamma*gamma/h/boltz/TK/sigma;
    I=beta_st*h/h/a/b*(4*h*sigma/3.14/h/a/b)^0.5*exp(-Delta_G/boltz/TK);
    Rp = 1.137*h*(I*v_st^2)^(1/3); % Calcite growth rate
end  
end



end


