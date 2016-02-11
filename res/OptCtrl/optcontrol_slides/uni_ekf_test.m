function f = uni_ekf_test
% EN530.603 Extended Kalman filtering of the unicycle with bearing and range measurements
% 
% M. Kobilarov , marin(at)jhu.edu

clear

%rng('default')
rng(10212)

S.bearing_only = 0;

% single beacon at (-2,2) : system is unobservable
%S.pbs = [-2;
%         2];    % beacon positions

% two beacons at (-2,2) and (2,2) : system is observable (two or more)
S.pbs = [-2, 2;
         2, 2];    % beacon positions

nb = size(S.pbs,2); % number of beacons

if S.bearing_only
  S.h = @b_h;    % bearing sensing
  S.r = nb;      % measurement dimension  
  S.R = .4*diag(repmat([.1], nb, 1));
else
  S.h = @br_h;   % bearing-reange sensing
  S.r = 2*nb;      % measurement dimension
  S.R = .4*diag(repmat([.1; .01], nb, 1));
end

S.n = 3;       % state dimension
S.f = @uni_f;  % mobile-robot dynamics

% timing
dt = .1;
%N = 2580;
N = 50;
T = dt*N;
S.dt = dt;

% noise models
S.Q = .3*dt*diag([.1 .1 .01]);


% initial mean and covariance
xt = [0; 0; pi/4]; % true state

P = .2*diag([1 1 .1]) % covariance
x = xt + sqrt(P)*randn(S.n, 1); % initial estimate with added noise

xts = zeros(S.n, N+1); % true states
xs = zeros(S.n, N+1);  % estimated states
Ps = zeros(S.n, S.n, N+1); % estimated covariances
ts = zeros(N+1,1); % times

zs = zeros(S.r, N);  % measurements

xts(:, 1) = xt;
xs(:, 1) = x;
Ps(:, :, 1) = P;
ts(1) = 0;

ds = zeros(S.n, N+1);  % errors
ds(:,1) = x - xt;

for k=1:N,
  u = dt*[2; 1];  % known controls
  
  xts(:,k+1) = S.f(xts(:,k), u, S) + sqrt(S.Q)*randn(3,1);  % true state

  [x,P] = ekf_predict(x, P, u, S);  % predict
  ts(k+1) = k*dt;
  
  z = S.h(xts(:,k+1), S) + sqrt(S.R)*randn(S.r,1);  % generate measurement  
  [x,P] = ekf_correct(x, P, z, S);  % correct
  
  xs(:,k+1) = x;
  Ps(:,:,k+1) = P;

  zs(:,k) = z;
  ds(:,k+1) = x - xts(:,k+1);  % actual estimate error
  ds(:,k+1) = fix_state(ds(:,k+1));
    
end

fig1=figure
fig2=figure

for k=1:N
  
  figure(fig1)
  plot(xts(1,1:k), xts(2,1:k), '--g','LineWidth',3)
  hold on
  plot(xs(1,1:k), xs(2,1:k), '-b','LineWidth',3)
  legend('true', 'estimated')

  % beacon
  plot(S.pbs(1,:), S.pbs(2,:), '*r');

  plotcov2(xs(1:2,k), 1.96^2*Ps(1:2,1:2,k));

%  quiver(xts(1,k), xts(2,k), .5*cos(xts(3,k)), .5*sin(xts(3,k)),'g');
%  quiver(xs(1,k), xs(2,k), .5*cos(xs(3,k)), .5*sin(xs(3,k)),'b');

  xlabel('x')
  ylabel('y')
  axis equal
  axis xy
  axis([-4 3 -4 4])

  saveas(gcf,['frames/frame' num2str(k,'%03d') '.jpg'])
  
  figure(fig2)
  subplot(2,1,1)
  hold off
  
  plot(ds(:,1:k)')

  xlabel('k')
  ylabel('meters or radians')
  legend('e_x','e_y','e_\theta')
  
  subplot(2,1,2)
  hold off
  plot(ts(1:k), reshape(sqrt(Ps(1,1,1:k)),k,1), ...
       ts(1:k), reshape(sqrt(Ps(2,2,1:k)),k,1), ...
       ts(1:k), reshape(sqrt(Ps(3,3,1:k)),k,1));
  legend('\sigma_x','\sigma_y','\sigma_\theta')
  
  xlabel('t')
  ylabel('meters or radians')
  
  drawnow
  
  if k==1
  
  end
  
%  saveas(frame,['frames/frame' num2str(k,'%03d') '.jpg'])

end

function [x, varargout] = uni_f(x, u, S)
% dynamical model of the unicycle
c = cos(x(3));
s = sin(x(3));

x = [x(1) + c*u(1);
     x(2) + s*u(1);
     x(3) + u(2)];

x = fix_state(x, S);

if nargout > 1
  % F-matrix
  varargout{1} = [1, 0, -s*u(1);
                  0, 1, c*u(1); 
                  0 0 1];
end


function [y, varargout] = br_h(x, S)

p = x(1:2);

y = [];
H = [];
for i=1:size(S.pbs, 2)
  pb = S.pbs(:, i); %i-th beacon
  d = pb - p;
  r = norm(d);
  
  th = fangle(atan2(d(2), d(1)) - x(3));
  y = [y; th; r];
  
  if nargout > 1
    % H-matrix
    H = [H;
         d(2)/r^2, -d(1)/r^2, -1;
         -d'/r, 0];
  end
end

if nargout > 1  
  varargout{1} = H;
end


function [y, varargout] = b_h(x, S)

p = x(1:2);

y = [];
H = [];
for i=1:size(S.pbs, 2)
  pb = S.pbs(:, i); %i-th beacon
  d = pb - p;
  r = norm(d);
  
  th = fangle(atan2(d(2), d(1)) - x(3));
  y = [y; th];
  
  if nargout > 1
    % H-matrix
    H = [H;
         d(2)/r^2, -d(1)/r^2, -1];
  end
end

if nargout > 1  
  varargout{1} = H;
end



function [x,P] = ekf_predict(x, P, u, S)

[x, F] = S.f(x, u, S);
x = fix_state(x, S);  % fix any [-pi,pi] issues
P = F*P*F' + S.Q;


function [x,P] = ekf_correct(x, P, z, S)

[y, H] = S.h(x, S);
P = P - P*H'*inv(H*P*H' + S.R)*H*P;
K = P*H'*inv(S.R);

e = z - y;
e = fix_meas(e, S);  % fix any [-pi,pi] issues
x = x + K*e;


function x = fix_state(x, S)
x(3) = fangle(x(3));

function z = fix_meas(z, S)
z
S.r
for i=1:size(S.pbs,2)
  if S.bearing_only
    z(i) = fangle(z(i));
  else
    z(2*i-1) = fangle(z(2*i-1));
  end
end

function a = fangle(a)
% make sure angle is between -pi and pi
a = mod(a,2*pi);
if a < -pi
  a = a + 2*pi;
else
  if a > pi
    a = a - 2*pi;
  end
end