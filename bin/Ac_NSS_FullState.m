%function [t,x,y_filter,u,delu,S] = Ac_NSS_FullState()
clear all;
close all;
tfinal=5;
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
[t,x]=ode45(@(t,x) NSS_overall_f( t,x, S), 0:0.001:tfinal, [S.x0;delx0;xf0;xf0;p0vec], options); 

[t_sim,x_sim]=ode45(@(t,x) sim_f( t,x, S), 0:0.001:tfinal, [S.x0-[1;0];S.x0-[1;0]], options_sim); 

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
y_sim_temp=lsim(sys,delu,t,S.x0-[1;0]);
y_sim_temp=lsim(sys,delu_sim,t,S.x0-[1;0]);

y_sim=lsim(H1,y_sim_temp,t);
%figure;
plot(t,y_sim,'g.',t,y_filter,'r-')
legend('Simulated signal','Actual Signal')
title('Simulates signal Vs Actual Signal')
xlabel('time') % x-axis label
ylabel('Output Signal after filtering') % y-axis label





