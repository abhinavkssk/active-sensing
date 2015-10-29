function [t,x,y_filter,u,delu,S] = Ac_SS_Fullstate()
close all;
tfinal=10;
S.k=1;

S.m=1;
S.b=0.1;
S.freq=5;
S.omega=2*pi*S.freq;
S.g=@sense_f;
S.ff=0.5;%Filter Factor



% initial state
x0 = [1;0];
S.x0 = x0;

xd= [3;0];
S.xd=xd;


options = odeset('RelTol',1e-4,'AbsTol',1e-4*ones(1,9));
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



[t,x]=ode45(@(t,x) overall_f( t,x, S), 0:0.001:30, [S.x0;S.x0-[1;0];zeros(5,1)], options); 

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
figure;


%plot(t,x(:,3),'g',t,delx_actual(:,1),'r');
%legend('delta x1 est','delta x2 actual')
%title('del est vs delx actual')
%xlabel('time') % x-axis label
%ylabel('X1_est X1_actual') % y-axis label


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


end




function dx = overall_f( t,x, S)
S.xstar=[cos(S.omega*t) -S.omega*sin(S.omega*t)];
S.ustar=-S.m*S.omega^2*cos(S.omega*t)-S.b*S.omega*sin(S.omega*t);
S.ystar=sense_nl(S.xstar,S);

x_sys=x(1:2);
delx_est=x(3:4);
x_fil=x(5:9);

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



dx_sys = A_sys*x_sys+B_sys*u;

ddelx_est=A_sys*delx_est+B_sys*delu+S.kfgain*(y_filter-S.C_sys*delx_est);

dx_fil=A_f*x_fil+B_f*dely_mod;

dx=[dx_sys;ddelx_est;dx_fil];
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

