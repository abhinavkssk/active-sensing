function f = int_kf_P1
%Farshid Alambeigi
%
clc
clear all
close all


%Timing Parameters
dt = 1;   % time-step
N = 10;   % total time-steps
T = N*dt;  % final time

%Noise Terms
S.n = (1.5*10^-5)^2; %measurement noise variance
S.u = (3*10^-9)^2;   %process noise of u variance
S.v = (3*10^-6)^2;  %process noise of v variance

%F matrix
S.F = [1 -dt;
       0 1];

%G matrix
S.G = [dt; 
       0];

%Q  matrix
S.Q = [S.v*dt+dt^3/3*S.u, -dt^2/2*S.u;
           -dt^2/2*S.u, dt*S.u];
       
%Noise models for true dynamics
S.Qn=diag([-sqrt(S.v) sqrt(S.u)]);

%R matrix
S.R = S.n;

%H matrix: theta measurment
S.H = [1, 0];

%Initial estimate of mean and covariance
x = [0; 1.7e-7]; 
P = diag([1e-4, 1e-12]); 

xts = zeros(2, N+1); % true states
xs = zeros(2, N+1);  % estimated states
Ps = zeros(2, 2, N+1); % estimated covariances

zs = zeros(1, N);  % estimated state

xts(:,1) = x;
xs(:,1) = x;
Ps(:,:,1) = P;

for k=1:N
  u = 0.02+sqrt(S.v)*randn(1)+xts(2,k); %Pick some known control
  xts(:,k+1) = S.F*xts(:,k) + S.G*u +S.Qn*randn(2,1);  %True state
  z = xts(1,k+1) + sqrt(S.n)*randn;   %Generate random measurement 
  
  [x,P] = kf_predict(x,P,u,S);  %Prediction
  [x,P] = kf_correct(x,P,z,S);  %Correction
  
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

% 96% confidence intervals of the estimated position
plot(xs(1,:) + 1.96*reshape(sqrt(Ps(1,1,:)),N+1,1)', '-k')
plot(xs(1,:) - 1.96*reshape(sqrt(Ps(1,1,:)),N+1,1)', '-k')


function [x,P] = kf_predict(x, P, u, S)

x = S.F*x + S.G*u;
P = S.F*P*S.F' + S.Q;

function [x,P] = kf_correct(x, P, z, S)

K = P*S.H'*inv(S.H*P*S.H' + S.R);
P = (eye(length(x)) - K*S.H)*P;
x = x + K*(z - S.H*x);
