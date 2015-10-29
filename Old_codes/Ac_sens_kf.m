function [t,x,y_filter,u] = Ac_sens_kf()
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

[t,x]=ode45(@(t,x) overall_f( t,x, S), [0 10], S.x0, options); 
u=t;
C=sense_lin(t,x,S);
xt=x';

for i=1:size(t,1)
    y_sys(i)=C(i,:)*xt(:,i);
    y_mod(i)=y_sys(i)*sin(S.omega*t(i));
end 
       
[ze,pe,ke] = ellip(5,3,30,S.omega/10,'s');
[be,ae] = zp2tf(ze,pe,ke);
H1=tf(be,ae);

tq = 0:0.001:10;
vq = interp1(t,y_mod,tq);
y_filter=lsim(H1,vq,tq);

H=-tf([3*S.omega],[2*S.m 2*S.b 0]);
y_sim=lsim(H*H1,u,t);

plot(t,y_sim,'.',0:0.001:10,y_filter,'-')

end

function dx = overall_f( t,x, S)

[~,A_sys,B_sys]=body_f(x,S);

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
