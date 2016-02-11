function f = uni_test1
% Extended Kalman filtering of the unicycle with bearing and range measurements

rng('default')

S.f = @uni_f;  % mobile-robot dynamics
S.h = @br_h;  % bearing-reange sensing
S.n = 4;      % state dimension
S.r = 2;      % measurement dimension

S.p0 = [0; 2];    % beacon position

% timing
dt = .1;
N = 50;
T = dt*N;
S.dt = dt;

% noise models
S.Q = .1*dt*dt*diag([.1 .1 .1,.001]);
S.R = .01*diag([.5 1]);%??

% initial mean and covariance
xt = [0; 0; 0; 1]; % true state

P = 10*.01*diag([2 2 2 5]) % covariance
x = xt + sqrt(P)*randn(S.n, 1); % initial estimate with added noise

xts = zeros(S.n, N+1); % true states
xs = zeros(S.n, N+1);  % estimated states
Ps = zeros(S.n, S.n, N+1); % estimated covariances

zs = zeros(S.r, N);  % measurements

xts(:, 1) = xt;
xs(:, 1) = x;
Ps(:, :, 1) = P;

ds = zeros(S.n, N+1);  % errors
ds(:,1) = x - xt;

for k=1:N,
  u = dt*[2; 1];  % known controls
  
  xts(:,k+1) = S.f(xts(:,k), u, S) + sqrt(S.Q)*randn(S.n,1);  % true state

  [x,P] = ekf_predict(x, P, u, S);  % predict
  
  z = S.h(xts(:,k+1), S) + sqrt(S.R)*randn(S.r,1);  % generate measurement
  z(1) = fangle(z(1));
  [x,P] = ekf_correct(x, P, z, S);  % correct
  
  xs(:,k+1) = x;
  Ps(:,:,k+1) = P;

  zs(:,k) = z;
  ds(:,k+1) = x - xts(:,k+1);  % actual estimate error
end

subplot(1, 2, 1)

plot(xts(1,:), xts(2,:), '--g','LineWidth',3)
hold on
plot(xs(1,:), xs(2,:), '-b','LineWidth',3)
legend('true', 'estimated')

xlabel('x')
ylabel('y')
axis equal
axis xy

% beacon
plot(S.p0(1), S.p0(2), '*r');

for k=1:5:N
  plotcov2(xs(1:2,k+1), Ps(1:2,1:2,k+1));
end

subplot(1,2,2)

plot(ds')

mean(sqrt(sum(ds.*ds, 1)))
xlabel('k')
ylabel('meters or radians')
legend('e_x','e_y','e_\theta','e_r')

figure, hold on;
plot(dt*(1:N+1),xs(4,:),'b');
plot(dt*(1:N+1),xts(4,:),'r');
xlabel('time(sec)');
ylabel('Radius(m)');
legend('Rest','Rtrue');

function [x, varargout] = uni_f(x, u, S)
% dynamical model of the unicycle modified
% x = [x,y,theta,radius]
% u = [sigma(commanded wheel vel), rotational angl vel];
c = cos(x(3));
s = sin(x(3));

x = [x(1) + c*x(4)*u(1);
     x(2) + s*x(4)*u(1);
     x(3) + u(2);
     x(4)];

if nargout > 1
  % F-matrix
  varargout{1} = [1, 0, -s*x(4)*u(1), c*u(1);
                  0, 1, c*x(4)*u(1), s*u(1); 
                  0, 0, 1, 0;
                  0, 0, 0, 1];  
end


function [y, varargout] = br_h(x, S)

p = x(1:2);
px = p(1);
py = p(2);

d = S.p0 - p;
r = norm(d);

th = fangle(atan2(d(2), d(1)) - x(3));

y = [th; r];

if nargout > 1
  % H-matrix
  varargout{1} = [d(2)/r^2, -d(1)/r^2, -1, 0;
                  -d'/r, 0, 0];
end


function [x,P] = ekf_predict(x, P, u, S)

[x, F] = S.f(x, u, S);
P = F*P*F' + S.Q;


function [x,P] = ekf_correct(x, P, z, S)

[y, H] = S.h(x, S);

K = P*H'*inv(H*P*H' + S.R);
P = (eye(S.n) - K*H)*P;

x = x + K*fangle(z-y);


function a = fangle(a)
% make sure angle is between -pi and pi
if a < -pi
  a = a + 2*pi
else
  if a > pi
    a = a - 2*pi
  end
end