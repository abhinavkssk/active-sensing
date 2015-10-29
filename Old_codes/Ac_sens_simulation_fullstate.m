function [t,x,S] = Ac_sens_simulation_fullstate()
close all;
S.tf=10;
S.k=1;

S.m=1;
S.b=0.1;
S.freq=5;
S.omega=2*pi*S.freq;
S.g=@sense_f;
S.ff=0.5;%Filter Factor

S.ycov = (1.5e-5)^2;


% initial state
x0 = [0;0];

S.x0 = x0;

X0=zeros(7,1);
X0(1)=S.x0(1);
X0(2)=S.x0(2);
X0(5)=0;
X0(7)= 0;

options = odeset('RelTol',1e-2,'AbsTol',1e-2*ones(1,7));

S.lqrPf = diag([200;200]);

S.lqrQ = diag([2; 2]);

S.lqrR = 10;

S.xd=[3;3];



sf = - S.lqrPf*S.xd;

csf = [P2c(S.lqrPf); sf];

[ts, css] = ode45(@(tf, csf) ricattiS(tf,csf,S), [S.tf:-0.01:0], csf);

css= flipdim(css,1);
S.lqrS11 = css(:, 1)';
S.lqrS12 = css(:, 2)';
S.lqrS22 = css(:, 3)';
S.lqrs1 = css(:, 4)';
S.lqrs2 = css(:, 5)';

[t,x]=ode45(@(t,x) overall_f( t,x, S), [0 S.tf], X0, options); 

plot(t,x(:,1:2))
figure;


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
