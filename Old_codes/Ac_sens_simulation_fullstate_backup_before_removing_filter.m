function [t,x,y_filter,u] = Ac_sens_simulation_fullstate()
close all;
tf=10;
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

X0=zeros(17,1);
X0(1)=S.x0(1);
X0(2)=S.x0(2);
options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4]);

S.lqrPf = diag([200;200]);

S.lqrQ = diag([2; 2]);

S.lqrR = 10;

S.xd=[3;3];



sf = - S.lqrPf*S.xd;

csf = [P2c(S.lqrPf); sf];

[ts, css] = ode45(@(tf, csf) ricattiS(tf,csf,S), [tf:-0.001:0], csf);

css= flipdim(css,1);
S.lqrS11 = css(:, 1)';
S.lqrS12 = css(:, 2)';
S.lqrS22 = css(:, 3)';
S.lqrs1 = css(:, 4)';
S.lqrs2 = css(:, 5)';

[t,x]=ode45(@(t,x) overall_f( t,x, S), [0 tf], X0, options); 



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
S.xstar=[cos(S.omega*t) -S.omega*sin(S.omega*t)];
S.ustar=-S.m*S.omega^2*cos(S.omega*t)-S.b*S.omega*sin(S.omega*t);

[~,A_sys,B_sys]=body_f(x,S);
C=sense_lin(t,x,S);

[ze,pe,ke] = ellip(5,3,30,S.omega/10,'s');
[be,ae] = zp2tf(ze,pe,ke);
[A_f,B_f,C_f,~] = tf2ss(be,ae);
A_del_filter=[A_sys zeros(2,5);sin(S.omega*t)*B_f*C A_f];
B_del_filter=[B_sys ;0*B_f];
C_del_filter=[0*C C_f];

x_sys=x(1:2);
x_fil=x(3:4);
delx_est=x(5:6);
x_fil_est=x(7:8);
delx_filter_est=[delx_est;x_fil_est];

p11=x(9);p12=x(10);p22=x(11);
p_mat=[ p11 p12;p12 p22];

dp_mat=A_del_filter*p_mat+p_mat*A_del_filter'-p_mat*C_del_filter'*inv(S.ycov)*C_del_filter*p_mat;
dp_vec=P2c(dp_mat);
dp11=dp_vec(1);dp12=dp_vec(2);dp22=dp_vec(3);
kal_K=p_mat*C_del_filter'*inv(S.ycov);



%x_est=delx_est+S.xstar;

tq = 0:0.001:tf;
s11 = interp1(t,S.lqrS11,tq);
s12 = interp1(t,S.lqrS12,tq);
s22 = interp1(t,S.lqrS22,tq);
s1 = interp1(t,S.lqrs1,tq);
s2 = interp1(t,S.lqrs2,tq);
P = c2P([s11; s12;s22]);
s = [s1;s2];
l = P*delx_est + s;
delu = -inv(S.lqrR)*B_sys'*l;

u=delu+S.ustar;

dx_sys = A_sys*x_sys+B_sys*u;

del_x=x_sys-S.xstar;
delx_filter=[del_x;x_fil];
y_del_filter=C_del_filter*delx_filter;

ddelx_filter_est=A_del_filter*delx_filter_est+B_del_filter*delu+kal_K*(y_del_filter-C_del_filter*delx_filter_est);

ddel_x_est=ddelx_filter_est(1:2);
dx_filter_est=ddelx_filter_est(3:4);

dx=[dx_sys;dx_filter;ddel_x_est;dx_filter_est;dp11;dp12;dp22];
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
