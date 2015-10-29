function [t,x,u] = Ac_sense_SteadyState( )
close all;
S.k=1;

S.m=1;
S.b=0.1;
S.freq=5;
S.omega=2*pi*S.freq;
S.g=@sense_f;
S.ff=0.5;%Filter Factor

tf2ss([3*S.omega],[2*S.m 2*S.b 0]);
% initial state
x0 = [0;0];
xd= [3;0];
S.x0 = x0;
S.xd=xd;
X0=[S.x0;S.x0];
Xd=[S.xd;S.xd];
options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4 1e-4]);

S.lqrN=[0;0];
S.lqrQ = diag([2; 2]);

S.lqrR = 10;
[~,S.A_sys,S.B_sys]=body_f(x0,S);
S.C_sys=[-47.1239 0];
S.D_sys=0;

[S.lqrgain,~,~] = lqr(S.A_sys,S.B_sys,S.lqrQ,S.lqrR,S.lqrN);

S.kfR = 1e-2;
S.kfQ=1e-2;
S.kfN=0;
sys=ss(S.A_sys,S.B_sys,S.C_sys,S.D_sys);
[~,S.kfgain,~] = kalman(sys,S.kfQ,S.kfR,S.kfN);

[t,x]=ode45(@(t,x) overall_f( t,x, S), 0:0.001:10, X0-Xd, options); 


for i=1:size(t,1)
u(i)=-S.lqrgain*x(i,3:4)';
end

x=x+repmat(Xd',size(t,1),1);
x=x(:,1:2);
plot(t,x,t,u);
xstar=[cos(S.omega*t) -S.omega*sin(S.omega*t)];
total_x=x+xstar;

for i=1:size(t,1)
ystar(i)=sense_nl(xstar(i,:),S);
total_y(i)=sense_nl(total_x(i,:),S);
end

y_sys_nl=total_y-ystar;
C=sense_lin(t,x,S);
xt=x';
for i=1:size(t,1)
y_sys(i)=C(i,:)*xt(:,i);
y_mod(i)=y_sys(i)*sin(S.omega*t(i));
y_mod_nl(i)=y_sys_nl(i)*sin(S.omega*t(i));
end
%figure;
%plot(t,y_sys_nl,'g.',t,y_sys,'b-')    
       


[ze,pe,ke] = ellip(5,3,30,S.omega/10,'s');
[be,ae] = zp2tf(ze,pe,ke);
H1=tf(be,ae);

%tq = 0:0.001:10;
%vq = interp1(t,y_mod,tq);
y_filter=lsim(H1,y_mod,t);
y_filter_nl=lsim(H1,y_mod_nl,t);

H=-tf([3*S.omega],[2*S.m 2*S.b 0])
y_sim_no_fildyn=lsim(H,u,t);
y_sim=lsim(H*H1,u,t);



figure;
plot(t,y_sim_no_fildyn,'g.',t,y_filter,'b-')

figure;
plot(t,y_sim,'g.',t,y_filter_nl,'r-')



end

function dx = overall_f( t,x, S)

A_sys=S.A_sys;
B_sys=S.B_sys;
C_sys=S.C_sys;
K=S.lqrgain;
L=S.kfgain;

dx = [A_sys -B_sys*K;
    L*C_sys  A_sys-B_sys*K-L*C_sys]*x;
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
gdash_val=3;
C=[-gdash_val*S.omega*sin(S.omega*t) g(cos(S.omega*t))];
end

function y = sense_nl(x,S)
g=S.g;

y=g(x(1))*x(2);


end

function y= sense_f(x)
y=3*x+5;
end

