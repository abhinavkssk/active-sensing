function [t,x,y_filter,u,delu,S] = Ac_NSS_FullState()
close all;
tfinal=10;
S.k=1;

S.m=1;
S.b=0.1;
S.freq=5;
S.omega=2*pi*S.freq;
S.g=@sense_f;
S.ff=0.5;%Filter Factor

S.ycov = 1.5e-4;
S.xcov=diag([1e-2 3e-2 1e-4 1e-4 1e-4 1e-4 1e-4]);
% initial state
x0 = [1;0];
S.x0 = x0;

xd= [3;0];
S.xd=xd;


options = odeset('RelTol',1e-4,'AbsTol',1e-4*ones(1,42));
options_sim = odeset('RelTol',1e-4,'AbsTol',1e-4*ones(1,4));




S.lqrN=[0;0];
S.lqrQ = diag([2; 2]);

S.lqrR = 10;
[~,S.A_sys,S.B_sys]=body_f(x0,S);
S.C_sys=[-47.1239 0];
S.D_sys=0;

[S.lqrgain,~,~] = lqr(S.A_sys,S.B_sys,S.lqrQ,S.lqrR,S.lqrN);

S.kfR = 1e-4;
S.kfQ=1e-4;
S.kfN=0;
sys=ss(S.A_sys,S.B_sys,S.C_sys,S.D_sys);
[~,S.kfgain,~] = kalman(sys,S.kfQ,S.kfR,S.kfN);


delx0=S.x0-[1;0];
xf0=zeros(5,1);
p0=[delx0;xf0;]*[delx0;xf0]';
p0vec=P2c(diag([1e-2 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]));
[t,x]=ode45(@(t,x) overall_f( t,x, S), 0:0.001:5, [S.x0;delx0;xf0;xf0;p0vec], options); 

[t_sim,x_sim]=ode45(@(t,x) sim_f( t,x, S), 0:0.001:30, [S.x0-[1;0];S.x0-[1;0]], options_sim); 

for i=1:size(t_sim,1)
delu_sim(i)=-S.lqrgain*(x_sim(i,3:4)'-S.xd);
end
for i=1:size(t,1)
delu(i)=-S.lqrgain*(x(i,3:4)'-S.xd);
end
S.ustar=-S.m*S.omega^2*cos(S.omega*t)-S.b*S.omega*sin(S.omega*t);
u=delu+S.ustar';
x=x(:,1:4);
%x=x+repmat([S.xd;S.xd]',size(t,1),1);
xstar=[cos(S.omega*t) -S.omega*sin(S.omega*t)];
x1star=cos(S.omega*t);

plot(t,x(:,3),t,x(:,1)-x1star);
legend('Est delx ','Actual delx')
title('Position x1 Vs t')
xlabel('time') % x-axis label
ylabel('X1') % y-axis label
%figure;


%plot(t,x(:,3),'g',t,delx_actual(:,1),'r');
%legend('delta x1 est','delta x2 actual')
%title('del est vs delx actual')
%xlabel('time') % x-axis label
%ylabel('X1_est X1_actual') % y-axis label

%{
for i=1:size(t,1)
ystar(i)=sense_nl(xstar(i,:),S);
total_y(i)=sense_nl(x(i,:),S);
end

y_sys_nl=total_y-ystar;
for i=1:size(t,1)
y_mod_nl(i)=y_sys_nl(i)*sin(S.omega*t(i));
end

[ze,pe,ke] = ellip(5,3,30,S.omega/10,'s');
[be,ae] = zp2tf(ze,pe,ke);
H1=tf(be,ae);

y_filter=lsim(H1,y_mod_nl,t);

%H=-tf([3*S.omega],[2*S.m 2*S.b 0]);
y_sim_temp=lsim(sys,delu,t);

y_sim=lsim(H1,y_sim_temp,t);
%figure;
plot(t,y_sim,'g.',t,y_filter,'r-')
legend('Simulated signal','Actual Signal')
title('Simulates signal Vs Actual Signal')
xlabel('time') % x-axis label
ylabel('Output Signal after filtering') % y-axis label
%}

end




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
function dx = sim_f( t,x, S)

A_sys=S.A_sys;
B_sys=S.B_sys;
C_sys=S.C_sys;
K=S.lqrgain;
L=S.kfgain;
x_sys=x(1:2);
x_est=x(3:4);
u=-K*(x_est-S.xd);
dx_sys = A_sys*x_sys +B_sys*u;
 dx_est= A_sys*x_est +B_sys*u+L*C_sys*(x_sys-x_est);%  A_sys-B_sys*K-L*C_sys]*x;
 
 dx=[dx_sys;dx_est];
end
function [x,A,B] = body_f( x, S)
b = S.b;
m=S.m;
A = [ 0 1;0 -b/m];

B = [0;1/m;];
end

function C = sense_lin(t,x,S)
g=S.g;

syms xsym;
gsym=@(xsym)g(xsym);

gd=diff(gsym(xsym));
gfun = matlabFunction(gd);
try
  gdash_val=gfun(x);
catch ME
  if (strcmp(ME.identifier ,'MATLAB:TooManyInputs'))     
   gdash_val=gfun();
  end
end

C=[-gdash_val*S.omega*sin(S.omega*t) g(cos(S.omega*t))];
end

function y = sense_nl(x,S)
g=S.g;

y=g(x(1))*x(2);


end

function y= sense_f(x)
y=3*x+5;
end

function c = P2c(P)


c = [P(1,1); P(1,2);P(1,3); P(1,4);P(1,5); P(1,6);P(1,7); P(2,2); P(2,3); P(2,4); P(2,5); P(2,6); P(2,7); P(3,3); P(3,4); P(3,5); P(3,6); P(3,7)...
    ; P(4,4); P(4,5); P(4,6); P(4,7); P(5,5); P(5,6); P(5,7); P(6,6); P(6,7); P(7,7);];
end

function P = c2P(c)
    c11=c(1);c12=c(2);c13=c(3);c14=c(4);c15=c(5);c16=c(6);c17=c(7);
    c21=c12;c22=c(8);c23=c(9);c24=c(10);c25=c(11);c26=c(12);c27=c(13);
    c31=c13;c32=c23;c33=c(14);c34=c(15);c35=c(16);c36=c(17);c37=c(18);
    c41=c14;c42=c24;c43=c34;c44=c(19);c45=c(20);c46=c(21);c47=c(22);
    c51=c15;c52=c25;c53=c35;c54=c45;c55=c(23);c56=c(24);c57=c(25);
    c61=c16;c62=c26;c63=c36;c64=c46;c65=c56;c66=(26);c67=c(27);
    c71=c17;c72=c27;c73=c37;c74=c47;c75=c57;c76=c67;c77=c(28);


P = [c11,c12,c13,c14,c15,c16,c17;
    c21,c22,c23,c24,c25,c26,c27;
    c31,c32,c33,c34,c35,c36,c37;
    c41,c42,c43,c44,c45,c46,c47;
    c51,c52,c53,c54,c55,c56,c57;
    c61,c62,c63,c64,c65,c66,c67;
    c71,c72,c73,c74,c75,c76,c77;];
end

