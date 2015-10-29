function [t,x,y_filter,u] = Ac_sens_simulation()
close all;
S.k=1;

S.m=1;
S.b=0.1;
S.freq=5;
S.omega=2*pi*S.freq;
S.g=@sense_f;
S.ff=0.5;%Filter Factor

% initial state
x0 = [0;0];

S.x0 = x0;

X0=[S.x0;0;0];

options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4]);

[t,x]=ode45(@(t,x) overall_f( t,x, S), 0:0.001:10, S.x0, options); 
u=t;

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
%plot(t,y_sys_nl,'g.',t,y_sys,'b-')    
       


[ze,pe,ke] = ellip(5,3,30,S.omega/10,'s');
[be,ae] = zp2tf(ze,pe,ke);
H1=tf(be,ae);

%tq = 0:0.001:10;
%vq = interp1(t,y_mod,tq);
y_filter=lsim(H1,y_mod,t);
y_filter_nl=lsim(H1,y_mod_nl,t);

H=-tf([3*S.omega],[2*S.m 2*S.b 0])
y_sim=lsim(H*H1,u,t);

figure;
plot(t,y_sim,'g.',t,y_filter,'b-')

figure;
plot(t,y_sim,'g.',t,y_filter_nl,'r-')
%{

figure;
plot(t,x(:,1),'-',t,x(:,2),'-.',t,y_filter,'.')
figure;
t_unif=(0:0.01:10)';
vq1 = sin(t_unif);
lsim(tf([2 5 1],[1 2 3]),vq1,t_unif);
%}
end

function dx = overall_f( t,x, S)

[~,A_sys,B_sys]=body_f(x,S);
%C_sys=sense_lin(t,x,S);
%[a,b,c,d]=butter(2,0.1,'low','s');

%A=[A_sys zeros(2);b*C_sys*sin(S.omega*t) a];
%B=[B_sys;zeros(2,1)];
%{
Controlmat=rank([B A*B A*A*B A*A*A*B])
if Controlmat<4
    pause
end
%}
u=t;
dx = A_sys*x+B_sys*u;
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
