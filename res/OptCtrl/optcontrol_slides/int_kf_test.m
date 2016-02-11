function f = int_kf_test
% Kalman filtering of the double integrator with position measurements

% timing
dt = .1;   % time-step
N = 100;   % total time-steps
T = N*dt;  % final time

% noise terms
S.q = .3;    % external disturbance variance
S.r = .9;    % measurement noise variance

% F matrix
S.F = [1 dt;
       0 1]

% G matrix
S.G = [dt^2/2; 
       dt];

% Q matrix
S.Q = S.q*[dt^3/3, dt^2/2;
           dt^2/2, dt];

% R matrix
S.R = S.r;

% H matrix
S.H = [1, 0];

% initial estimate of mean and covariance
x = [0; 0];
P = 10*diag([1 1]);

xts = zeros(2, N+1); % true states
xs = zeros(2, N+1);  % estimated states
Ps = zeros(2, 2, N+1); % estimated covariances

zs = zeros(1, N);  % estimated state

pms = zeros(1, N); % measured position

xts(:,1) = x;
xs(:,1) = x;
Ps(:,:,1) = P;

for k=1:N
  u = cos(k/N); % pick some known control

  xts(:,k+1) = S.F*xts(:,k) + S.G*(u + sqrt(S.q)*randn);  % true state

  [x,P] = kf_predict(x,P,u,S);  % prediction
  
  z = xts(1,k+1) + sqrt(S.r)*randn;   % generate random measurement 
  
  [x,P] = kf_correct(x,P,z,S);  % correction
  
  % record result
  xs(:,k+1) = x;
  Ps(:,:,k+1) = P;
  zs(:,k) = z;
end

plot(xts(1,:), '--', 'LineWidth',2)
hold on
plot(xs(1,:), 'g', 'LineWidth',2)
plot(2:N+1,zs(1,:), 'r', 'LineWidth',2)

legend('true', 'estimated','measured')

% 95% confidence intervals of the estimated position
plot(xs(1,:) + 1.96*reshape(sqrt(Ps(1,1,:)),N+1,1)', '-g')
plot(xs(1,:) - 1.96*reshape(sqrt(Ps(1,1,:)),N+1,1)', '-g')


function [x,P] = kf_predict(x, P, u, S)

x = S.F*x + S.G*u;
P = S.F*P*S.F' + S.Q;

function [x,P] = kf_correct(x, P, z, S)

K = P*S.H'*inv(S.H*P*S.H' + S.R);
P = (eye(length(x)) - K*S.H)*P;
x = x + K*(z - S.H*x);
