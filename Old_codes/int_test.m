function f = int_test
% Kalman filtering of the double integrator with position measurements

% timing
dt = 1;   % time-step
N = 30;   % total time-steps
T = N*dt;  % final time

% noise terms
S.qu = (3e-9)^2;    % external disturbance variance thetadot
S.qv = (3e-6)^2;    % external disturbance variance bias
S.qn = (1.5e-5)^2;    % measurement noise variance

% PHI matrix 
S.Phi = [1 -dt;
       0 1];

% G matrix
S.G = [dt; 
       0];

% Q matrix
S.Q = [S.qv*dt + S.qu*(dt^3/3), -S.qu*(dt^2/2);
           -S.qu*(dt^2/2), S.qu*dt];

% R matrix
S.R = S.qn;

% H matrix
S.H = [1, 0];

% initial estimate of mean and covariance
x = [0; 1.7e-7];
P = diag([1e-2; 1e-12]);
thetadot = 0.02; % given trajectory for true state theta 

xts = zeros(2, N+1); % true states
xs = zeros(2, N+1);  % estimated states
Ps = zeros(2, 2, N+1); % estimated covariances
nm = zeros(N+1,1);

zs = zeros(1, N);  % estimated state

pms = zeros(1, N); % measured position

xts(:,1) = x;
xs(:,1) = x;
Ps(:,:,1) = P;
nm(1) = norm(P);

for k=1:N
  xts(:,k+1) = xts(:,k) + [thetadot*dt;sqrt(S.qu)*randn];
  
  %generate u based on true state
  u = thetadot + xts(2,k+1) + sqrt(S.qv)*randn;

  [x,P] = kf_predict(x,P,u,S);  % prediction
  
  z = xts(1,k+1) + sqrt(S.qn)*randn;   % generate random measurement 
  
  [x,P] = kf_correct(x,P,z,S);  % correction
  
  % record result
  xs(:,k+1) = x;
  Ps(:,:,k+1) = P;
  zs(:,k) = z;
  nm(k+1) = norm(P);
end

plot(xts(1,:), 'x--', 'LineWidth',2)
hold on
plot(xs(1,:), 'gx-', 'LineWidth',2)
plot(dt*(2:N+1),zs(1,:), 'ro-', 'LineWidth',2)

xlabel('time(sec)');
ylabel('\theta(rad)');
legend('true', 'estimated','measured')

% 95% confidence intervals of the estimated position
plot(xs(1,:) + 1.96*reshape(sqrt(Ps(1,1,:)),N+1,1)', '-g')
plot(xs(1,:) - 1.96*reshape(sqrt(Ps(1,1,:)),N+1,1)', '-g')

figure;
plot(dt*(2:N+1),nm(2:N+1),'-');
legend('norm covariance');
xlabel('time(sec)');

figure, hold on;
error = xs - xts;
plot(dt*(1:N+1), error(1,:),'r');
plot(dt*(1:N+1), error(2,:),'b');
ylabel('rad or rad/s');
xlabel('time(s)');
legend('etheta','ebias');

function [x,P] = kf_predict(x, P, u, S)

x = S.Phi*x + S.G*u;
P = S.Phi*P*S.Phi' + S.Q;

function [x,P] = kf_correct(x, P, z, S)

K = P*S.H'*inv(S.H*P*S.H' + S.R);
P = (eye(length(x)) - K*S.H)*P;
x = x + K*(z - S.H*x);
