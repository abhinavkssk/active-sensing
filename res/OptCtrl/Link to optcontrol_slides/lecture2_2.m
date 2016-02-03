function res = lecture2_2()
% EN530.603 Applied Optimal Control
% Lecture#2, 2013/09/09
%
% Unconstrained nonlinear optimization basics (2 of 2)
% 
% Marin Kobilarov, marin(at)jhu.edu


% nonlinear function of two variables
syms x1 x2 real
x = [x1; x2]; % column vector of symbolic variables
%f = log(1 + 3*(x2 - (x1^3 - x1))^2 + (x1 - 4/3)^2);
f = x1*exp(-(x1^2 + x2^2)) + (x1^2 + x2^2)/20;

figure
ezsurfc(f,[-2 2])

figure
ezcontour(f,[-2 2],100)
hold on

% compute gradient and hessian symbolically
gradf = jacobian(f,x).' % column gradf
hessf = jacobian(gradf,x)

% handles for the function, its gradient, and Hessian
fx = matlabFunction(f,'vars',{x});
gx = matlabFunction(gradf,'vars',{x});
Hx = matlabFunction(hessf,'vars',{x});
fgHx = matlabFunction(f,gradf,hessf,'vars',{x});

reg = 1;  % whether to use regularized-Newton or not
x0 = [1.6; 1.9]; % reg-Newton and gradient converge; Newton is stuck
%x0=[-.4; 0.6]; % reg-Newton and gradient converge
% x0=[1; 0.1];  % gradient get stuck


% Test gradient method with Armijo rule
x=x0;
xs = x;

for k=1:100,
  x = gstep(fx, gx, x);  
  xs = [xs, x];
end

plot(xs(1,:), xs(2,:), '-ob')

disp(['Gradient method: x=[' num2str(x') '] f=' num2str(fx(x))]);

% Test Newton method with Armijo rule
x = x0;
xs = x;
for k=1:20,
  x = Nstep(fx, gx, Hx, x, reg);  
  xs = [xs, x];
end
plot(xs(1,:), xs(2,:), '-dm')
disp(['Newton method: x=[' num2str(x') '] f=' num2str(fx(x))]);


% run a standard algorithm
opts = optimoptions('fminunc','GradObj','on','Hessian','on');
[x, fval, exitflag, output] = fminunc(fgHx,x0,opts)
plot(x(1), x(2), '*g')
disp(['Fminunc: x=[' num2str(x') '] f=' num2str(fval)]);


function x = gstep(fx, gx, x)
% gradient step

% gradient
g = gx(x);

% test for convergence
if norm(g) < 1e-8
  return
end

% direction
d = -g;

% step-size
a = Armijo(fx, x, g, d);

% update state
x = x + a*d;


function x = Nstep(fx, gx, Hx, x, reg)
% Newton step

g = gx(x); % gradient
H = Hx(x); % Hessian

% test for convergence
if norm(g) < 1e-8
  return
end

% if regularized version
if reg
  % check for positive definiteness  
  e = eig(H);
  if (min(e) < 0)
    % apply a trust-region like Hessian modification
    H = H + eye(2)*(abs(min(e)) + .001);
  end
end

% direction, i.e. d=-inv(H)*g;
d = -H\g;

%step-size
a = Armijo(fx, x, g, d);  

% new state
x = x + a*d;    



function a = Armijo(fx, x, g, d)

sigma = .01;
beta = .25;
s = 1;

f = fx(x);
for m=0:1000
  a = beta^m*s;
  if (f - fx(x + a*d) > -sigma*beta^m*s*g'*d)
    return
  end
end
