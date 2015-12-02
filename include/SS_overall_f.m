

function dx = SS_overall_f( t,x, S)
t;
S.xstar=[cos(S.omega*t) -S.omega*sin(S.omega*t)];
S.ustar=-S.m*S.omega^2*cos(S.omega*t)-S.b*S.omega*sin(S.omega*t);
S.ystar=sense_nl(S.xstar,S);

x_sys=x(1:2);
delx_est=x(3:4);
x_fil=x(5:9);
x_fil_est=x(10:14);

A_sys=S.A_sys;
B_sys=S.B_sys;


[ze,pe,ke] = ellip(5,3,30,S.omega/10,'s');
[be,ae] = zp2tf(ze,pe,ke);
[A_f,B_f,C_f,~] = tf2ss(be,ae);
y_filter=C_f*x_fil;




delu = -S.lqrgain*(delx_est-S.xd);

u=delu+S.ustar;

y=sense_nl(x_sys,S);


dely=y-S.ystar;
dely_mod=dely*sin(S.omega*t);

delx_fil_est=[delx_est;x_fil_est];

dx_sys = A_sys*x_sys+B_sys*u;

ddelx_fil_est=S.A_del_filter*delx_fil_est+S.B_del_filter*delu+S.kfgain*(y_filter-S.C_del_filter*delx_fil_est);

dx_fil=A_f*x_fil+B_f*dely_mod;

dx=[dx_sys;ddelx_fil_est(1:2);dx_fil;ddelx_fil_est(3:7)];
end











