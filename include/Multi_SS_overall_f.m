

function dx = Multi_SS_overall_f( t,x, S)
t;
S.xstar=[cos(S.omega*t) -S.omega*sin(S.omega*t)];
S.ustar=-S.m*S.omega^2*cos(S.omega*t)-S.b*S.omega*sin(S.omega*t);
S.ystar=sense_nl(S.xstar,S);

x_sys=x(1:2);
delx_est=x(3:4);
x_1fil=x(5:9);
x_1fil_est=x(10:14);
x_2fil=x(15:19);
x_2fil_est=x(20:24);
x_3fil=x(25:29);
x_3fil_est=x(30:34);

A_sys=S.A_sys;
B_sys=S.B_sys;


[ze,pe,ke] = ellip(5,3,30,S.omega/10,'s');
[be,ae] = zp2tf(ze,pe,ke);
[A_f,B_f,C_f,~] = tf2ss(be,ae);
y_filter=C_f*x_1fil;

y_3filter=repmat(y_filter,3,1);


delu = -S.lqrgain*(delx_est-S.xd);

u=delu+S.ustar;

y=sense_nl(x_sys,S);


dely=y-S.ystar;
dely_mod=dely*sin(S.omega*t);

delx_3fil_est=[delx_est;x_1fil_est;x_2fil_est;x_3fil_est];

dx_sys = A_sys*x_sys+B_sys*u;
%{
S.A_del_3filter*delx_3fil_est
S.B_del_3filter*delu
S.kfgain
y_3filter
S.C_del_3filter*delx_3fil_est
S.kfgain*(y_3filter-S.C_del_3filter*delx_3fil_est)
%}
ddelx_3fil_est=S.A_del_3filter*delx_3fil_est...
    +S.B_del_3filter*delu...
    +S.kfgain*(y_3filter-S.C_del_3filter*delx_3fil_est);

dx_1fil=A_f*x_1fil+B_f*dely_mod;
dx_2fil=A_f*x_2fil+B_f*dely_mod;
dx_3fil=A_f*x_3fil+B_f*dely_mod;

dx=[dx_sys;ddelx_3fil_est(1:2);dx_1fil;ddelx_3fil_est(3:7);dx_2fil;ddelx_3fil_est(8:12);dx_3fil;ddelx_3fil_est(13:17)];
end











