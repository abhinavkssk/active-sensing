function [dus, V, Vn, dV, a] = ddp(x0, us, S)
% Second-order numerical optimal control. The code computes
% the optimal control adjustment for a given dynamical system
%
% params:
% x0 - initial state
% us - m-N matrix with discrete controls
% S - problem data:
%     S.L  : handle to the cost function
%     S.f  : handle to the discrete dynamics functions
%     S.mu : regularizing constant (default is 0)
%     S.a  : initial step size (default is 1)
%     S.diff : difference function (default is minus, i.e. vector space)
%
% return:
%   dus: m-N matrix containing computed optimal change in control
%   V: current value function
%   Vn: new value function
%   dV: predicted change in value function
%   a: computed step-size along control search direction
%
%
% Note: this implementation is most closely related to second-order 
% metehods known as stage-wise Newton (SN) - Bertsekas, 2003 
% and differential-dynamic-programming (DDP), Mayne, 1966. 
% In this implementation second-order terms in the dynamics 
% are ignored which corresponds to the linear-quadratic-subproblem
% (LQS) approach (see also iterative-LQR (Todorov et al)).
%
% Disclaimer: the code is for education purposes only
%
% Author: Marin Kobilarov marin(at)jhu.edu

if ~isfield(S, 'diff')
  S.diff = @diff_def;
end

if ~isfield(S, 'mu')
  S.mu = 1e-4;
end

if ~isfield(S, 'a')
  S.a = 1;
end

if ~isfield(S, 'n')
  S.n = length(x0);
end

n = S.n;
m = size(us, 1);

N = size(us, 2);

Ps = zeros(n,n,N+1);
vs = zeros(n,N+1);

cs = zeros(m,N);
Ds = zeros(m,n,N);

dus = zeros(size(us));

% integrate trajectory and get terminal cost
xs = ddp_traj(x0, us, S);
[L, Lx, Lxx, Lu, Luu] = S.L(N+1, xs(:,end), [], S);


% initialize
V = L;
v = Lx;
P = Lxx;

dV = [0; 0];

Ps(:,:,N+1) = P;
vs(:,N+1) = v;

for k=N:-1:1,
  
  x = xs(:,k);
  u = us(:,k);
  
  [xn, A, B] = S.f(k, x, u, S);
  
  if isempty(A) || isempty(B)
    [A, B] = fd(S.f, k, x, u, S, 1e-6);    
  end  
  
  [L, Lx, Lxx, Lu, Luu] = S.L(k, x, u, S);
  
  V = V + L;
  
  Qx = Lx + A'*v;
  Qu = Lu + B'*v;
  Qxx = Lxx + A'*P*A;
  Quu = Luu + B'*P*B;
  Qux = B'*P*A;
  
  
  for count = 1:10
      Quum = Quu + S.mu*eye(m);
      [F, d] = chol(Quum);
      if d~= 0
        disp('[W] ddp: cholesky failed. Increasing mu:')
        S.mu = 10*S.mu;
      else
          break;
      end 
  end
  if d~=0
      disp('[W] ddp: cholesky failed. for mu:')
      disp(S.mu);
      dus = [];
      return 
  end    
  % control law is du = c + D*dx
  cD = -F\(F'\[Qu, Qux]);
  c = cD(:, 1);
  D = cD(:, 2:end);
  
  v = Qx + D'*Qu;
  P = Qxx + D'*Qux + D'*Quu*D;
  
  dV = dV + [c'*Qu; c'*Quu*c/2];
  
  vs(:, k) = v;
  Ps(:, :, k) = P;

  cs(:, k) = c; 
  Ds(:, :, k) = D; 

end

dV

s1 = .1;
s2 = .5;
%b1 = .25;
b1 = 0.25;
b2 = 2;

a = S.a;

% measured change in V
dVm = eps;

while dVm > 0

  % variation
  dx = zeros(n, 1);
  
  % varied x
  xn = x0;

  % new measured cost
  Vn = 0;% Should put J instead of V its kind of confusing
  
  for k=1:N,
    
    u = us(:,k);
    
    c = cs(:,k);
    D = Ds(:,:,k);

    du = a*c + D*dx;
    un = u + du;
    
    [Ln, Lx, Lxx, Lu, Luu] = S.L(k, xn, un, S);
    
    [xn, A, B] = S.f(k, xn, un, S);
    
    dx = S.diff(xs(:,k+1), xn);

    Vn = Vn + Ln;
    
    dus(:,k) = du;
  end
  
  [L, Lx, Lxx, Lu, Luu] = S.L(N+1, xn, [], S);
  Vn = Vn + L;
  
  dVm = Vn - V
  
  
%   if abs(dVm) < 1e-12
%       break
%   end
  if dVm > 0
    a = b1*a;
    disp(['ddp: decreasing a=' num2str(a)])
    continue    
  end
      
  dVp = [a; a*a]'*dV;
  
  r = dVm/dVp
  
  if r < s1
    a = max(b1*a,1e-16);
  else
    if r >= s2 
      a = b2*a;
    end
  end
  disp(['ddp: decreasing a=' num2str(a)])
  
end


function dx = diff_def(x, xn)
% default state difference 

dx = xn - x;


function [A, B] = fd(func, k, x, u, S, e)
% compute numerically the jacobians A=fx, B=fu of a given function f(k,x,u,S)

f = func(k, x, u, S);

n = length(x);
m = length(u);

En = eye(n);
Em = eye(m);

A = zeros(n, n);
B = zeros(n, m);

for j=1:n,
  A(:,j) = (func(k, x + e*En(:,j), u, S) - f)/e;
end

for j=1:m,
  B(:,j) = (func(k, x, u + e*Em(:,j), S) - f)/e;
end
