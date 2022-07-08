
function [Pa,Pm,Pb]=solveP(ka,kb,km,kbm,va,vb,vm,vbm,pa,pm,a,b,m,lmd)

dP=1e-8;
alpha=0.2;

for j=1:50
    
    P_A1=pa;
    P_M1=pm;
    P_B1=(1-P_A1-P_M1)/(1+lmd);
    u_AB1=ka*a*P_B1-va*P_A1;
    u_MB1=km*m*P_B1-vm*P_M1;
    rx1=u_MB1/u_AB1;
    x1=rx1/(1+rx1);
    u_BA1=kb*b*P_A1-vb*P_B1*(1-x1);
    u_BM1=kbm*b*P_M1-vbm*P_B1*x1;
    EA1=u_AB1-u_BA1;
    EM1=u_MB1-u_BM1;
    
    P_A2=pa+dP;
    P_M2=pm;
    P_B2=(1-P_A2-P_M2)/(1+lmd);
    u_AB2=ka*a*P_B2-va*P_A2;
    u_MB2=km*m*P_B2-vm*P_M2;
    rx2=u_MB2/u_AB2;
    x2=rx2/(1+rx2);
    u_BA2=kb*b*P_A2-vb*P_B2*(1-x2);
    u_BM2=kbm*b*P_M2-vbm*P_B2*x2;
    EA2=u_AB2-u_BA2;
    EM2=u_MB2-u_BM2;
    
    P_A3=pa;
    P_M3=pm+dP;
    P_B3=(1-P_A3-P_M3)/(1+lmd);
    u_AB3=ka*a*P_B3-va*P_A3;
    u_MB3=km*m*P_B3-vm*P_M3;
    rx3=u_MB3/u_AB3;
    x3=rx3/(1+rx3);
    u_BA3=kb*b*P_A3-vb*P_B3*(1-x3);
    u_BM3=kbm*b*P_M3-vbm*P_B3*x3;
    EA3=u_AB3-u_BA3;
    EM3=u_MB3-u_BM3;
    
    
    dEA_dPA=(EA2-EA1)/dP;
    dEM_dPA=(EM2-EM1)/dP;
    dEA_dPM=(EA3-EA1)/dP;
    dEM_dPM=(EM3-EM1)/dP;
    
    J=[dEM_dPA dEM_dPM;dEA_dPA dEA_dPM];
    
    dP_matrix=J\[EM1;EA1];
    pa=pa-alpha*dP_matrix(1);
    pm=pm-alpha*dP_matrix(2);
    
end

Pa=pa;
Pm=pm;
Pb=(1-pa-pm)/(1+lmd);

end