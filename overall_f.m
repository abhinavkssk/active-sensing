

function dx = overall_f( t,x, S)
t
S.xstar=[cos(S.omega*t) -S.omega*sin(S.omega*t)];
S.ustar=-S.m*S.omega^2*cos(S.omega*t)-S.b*S.omega*sin(S.omega*t);
S.ystar=sense_nl(S.xstar,S);

x_sys=x(1:2);
delx_est=x(3:4);
x_fil_est=x(5:9);
x_fil=x(10:14);
p_vec=x(15:42);

A_sys=S.A_sys;
B_sys=S.B_sys;



[ze,pe,ke] = ellip(5,3,30,S.omega/10,'s');
[be,ae] = zp2tf(ze,pe,ke);
[A_f,B_f,C_f,~] = tf2ss(be,ae);
y_filter=C_f*x_fil;

Ct=sense_lin(t,x,S);
A_del_filter=[A_sys zeros(2,5);sin(S.omega*t)*B_f*Ct A_f];
B_del_filter=[B_sys ;0*B_f];
C_del_filter=[0*Ct C_f];

p_mat=c2P(p_vec);
dp_mat=A_del_filter*p_mat+p_mat*A_del_filter'-p_mat*C_del_filter'*inv(S.ycov)*C_del_filter*p_mat+S.xcov;
dp_vec=P2c(dp_mat);

overall_kfgain=p_mat*C_del_filter'*inv(S.ycov);



delu = -S.lqrgain*(delx_est-S.xd);

u=delu+S.ustar;

y=sense_nl(x_sys,S);


dely=y-S.ystar;
dely_mod=dely*sin(S.omega*t);



dx_sys = A_sys*x_sys+B_sys*u;

overall_est=[delx_est;x_fil_est];
doverall_est=A_del_filter*overall_est+B_del_filter*delu+overall_kfgain*(y_filter-C_del_filter*overall_est);

dx_fil=A_f*x_fil+B_f*dely_mod;

dx=[dx_sys;doverall_est;dx_fil;dp_vec]
end











