function f=test(t,x,S)
[~,A_sys,B_sys]=body_f(x,S);

tq = 0:0.01:S.tf;
for i=1:size(t,1)
s11 = interp1(tq,S.lqrS11,t(i));
s12 = interp1(tq,S.lqrS12,t(i));
s22 = interp1(tq,S.lqrS22,t(i));
s1 = interp1(tq,S.lqrs1,t(i));
s2 = interp1(tq,S.lqrs2,t(i));
P = c2P([s11; s12;s22]);
s = [s1;s2];
x_req=[x(i,1);x(i,2)];
l = P*(x_req) + s;
u(i) = -inv(S.lqrR)*B_sys'*l;
end

x=x(:,1:2);
xstar=[cos(S.omega*t) -S.omega*sin(S.omega*t)];
for i=1:size(t,1)
ystar(i)=sense_nl(xstar(i,:),S);
total_y(i)=sense_nl(x(i,:),S);
end
y_sys_nl=total_y-ystar;

delx=x-xstar;

xt=delx';
C=sense_lin(t,delx,S);

for i=1:size(t,1)
y_sys(i)=C(i,:)*xt(:,i);
y_mod(i)=y_sys(i)*sin(S.omega*t(i));
y_mod_nl(i)=y_sys_nl(i)*sin(S.omega*t(i));
end
plot(t,y_sys_nl,'g.',t,y_sys,'b-')    
       
pause

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
end













function dcs = ricattiS(t,cs,S)

[~,A,B]=body_f(cs,S);
P = c2P(cs(1:3));
s = cs(4:5);
Q=S.lqrQ;
R=S.lqrR;
xd=S.xd;
dP = -P*A - A'*P - Q + P*B*inv(R)*B'*P;
ds = - (A' - P*B*inv(R)*B')*s + Q*xd;
dcs = [P2c(dP); ds];
end


function dx = overall_f( t,x, S)
S.xstar=[cos(S.omega*t); -S.omega*sin(S.omega*t)];
S.ustar=-S.m*S.omega^2*cos(S.omega*t)-S.b*S.omega*sin(S.omega*t);

[~,A_sys,B_sys]=body_f(x,S);
C=sense_lin(t,x,S);
C_af=[3*S.omega/2 0];

x_sys=x(1:2);
delx_est=x(3:4);

p11=x(5);p12=x(6);p22=x(7);
p_mat=[ p11 p12;p12 p22];

dp_mat=A_sys*p_mat+p_mat*A_sys'-p_mat*C_af'*inv(S.ycov)*C_af*p_mat;
dp_vec=P2c(dp_mat);
dp11=dp_vec(1);dp12=dp_vec(2);dp22=dp_vec(3);
kal_K=p_mat*C_af'*inv(S.ycov);



%x_est=delx_est+S.xstar;

tq = 0:0.01:S.tf;
s11 = interp1(tq,S.lqrS11,t);
s12 = interp1(tq,S.lqrS12,t);
s22 = interp1(tq,S.lqrS22,t);
s1 = interp1(tq,S.lqrs1,t);
s2 = interp1(tq,S.lqrs2,t);
P = c2P([s11; s12;s22]);
s = [s1;s2];
l = P*(S.xstar+delx_est) + s;
u = -inv(S.lqrR)*B_sys'*l;
delu=u-S.ustar;

dx_sys = A_sys*x_sys+B_sys*u;

del_x=x_sys-S.xstar; 

y_af=C_af*del_x;

ddel_x_est=A_sys*delx_est+B_sys*delu+kal_K*(y_af-C_af*delx_est);


dx=[dx_sys;ddel_x_est;dp11;dp12;dp22]
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
c = [P(1,1); P(1,2); P(2,2)];
end

function P = c2P(c)
P = [c(1), c(2);
     c(2), c(3)];
end
