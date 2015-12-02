function dx = Multi_sim_f( t,x, S)

A_sys=S.A_sys;
B_sys=S.B_sys;
C_sys=S.C_sys;
K=S.lqrgain;
L=S.kfgain_sim;
x_sys=x(1:2);
x_est=x(3:4);
u=-K*(x_est-S.xd);
dx_sys = A_sys*x_sys +B_sys*u;
 dx_est= A_sys*x_est +B_sys*u+L*C_sys*(x_sys-x_est);%  A_sys-B_sys*K-L*C_sys]*x;
 
 dx=[dx_sys;dx_est];
end