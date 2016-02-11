function f = uni_ekf_P3
% Farshid Alambeigi

clc
clear all
close all
clear

%rng('default')
rng(10212)  %Control the random number generator

S.bearing_only = 0;


%Two beacons at (-2,2) and (2,2) : system is observable (two or more)
S.pbs = [0, 0;
         5, 2];    %Beacon positions

nb = size(S.pbs,2); %Number of beacons

if S.bearing_only
  S.h = @b_h;    %Bearing sensing
  S.r = nb;      %Measurement dimension  
  S.R = .4*diag(repmat([.1], nb, 1));
else
  S.h = @br_h;   %Bearing-range sensing
  S.r = 2*nb;      %Measurement dimension
  S.R = .4*diag(repmat([.1; .01], nb, 1));
end

%S.h = @br_h; %Bearing-Range Sensing
%S.r = 2; %Measurement Dimension
S.n = 4;       %State dimension
S.f = @uni_f;  %Mobile-robot dynamics


%Timing Parameters
dt = .1;
N = 60;
T = dt*N;
S.dt = dt;

%Noise models
S.Q = dt^2*diag([.01 .01 .01 .0001]);    %Dynamics noise covariance


%Initial mean and covariance
xt = [0; 0; 0; 1];  %True state

P = diag([.01 .01 .01 .04]);   %Covariance
x = xt + sqrt(P)*randn(S.n, 1);     %Initial estimate with added noise

xts = zeros(S.n, N+1);  %True states
xs = zeros(S.n, N+1);   %Estimated states
Ps = zeros(S.n, S.n, N+1);  %Estimated covariances
ts = zeros(N+1,1);  %Times

zs = zeros(S.r, N);    %Measurements

xts(:, 1) = xt;
xs(:, 1) = x;
Ps(:, :, 1) = P;
ts(1) = 0;

ds = zeros(S.n, N+1);  %Errors
ds(:,1) = x - xt;


for k=1:N,
  u = dt*[2; 1];  %Known controls
  
  xts(:,k+1) = S.f(xts(:,k), u, S) + sqrt(S.Q)*randn(4,1);      %True state

  [x,P] = ekf_predict(x, P, u, S);      %Predict
  ts(k+1) = k*dt;
   z = S.h(xts(:,k+1), S) + sqrt(S.R)*randn(S.r,1);  %Generate measurement
    
  [x,P] = ekf_correct(x, P, z, S);  %Correct
  
  xs(:,k+1) = x;
  Ps(:,:,k+1) = P;

  zs(:,k) = z;
  ds(:,k+1) = x - xts(:,k+1);  %Actual estimate error
  ds(:,k+1) = fix_state(ds(:,k+1));
  dssize=size(ds,1)
end

subplot(2, 2, 1)
plot(xts(1,:), xts(2,:), '--r','LineWidth',3)
hold on
plot(xs(1,:), xs(2,:), '-.b','LineWidth',3)
legend('true', 'estimated')


xlabel('x')
ylabel('y')
axis equal
axis xy

%Beacon
plot(S.pbs(1,:), S.pbs(2,:), '*r');
grid on

for k=1:1:N
  plotcov2(xs(1:2,k), 1.96^2*Ps(1:2,1:2,k));
% if rem(k,2)==0
% quiver(xts(1,k), xts(2,k), .5*cos(xts(3,k)), .5*sin(xts(3,k)),'k');
% quiver(xs(1,k), xs(2,k), .5*cos(xs(3,k)), .5*sin(xs(3,k)),'g');
% end
end
quiver(xts(1,:), xts(2,:), .5*cos(xts(3,:)), .5*sin(xts(3,:)),'k');
quiver(xs(1,:), xs(2,:), .5*cos(xs(3,:)), .5*sin(xs(3,:)),'g');

subplot(2,2,2)
plot(ds','LineWidth',2)
grid on

mean(sqrt(sum(ds.*ds, 1)))
xlabel('k')
ylabel('meters or radians')
legend('e_x','e_y','e_\theta','e_r')

subplot(2,2,3)

plot(ts, reshape(sqrt(Ps(1,1,:)),N+1,1), ...
     ts, reshape(sqrt(Ps(2,2,:)),N+1,1), ...
     ts, reshape(sqrt(Ps(3,3,:)),N+1,1),...
     ts, reshape(sqrt(Ps(4,4,:)),N+1,1),'LineWidth',2);
 grid on
legend('\sigma_x','\sigma_y','\sigma_\theta','\sigma_r')

xlabel('t')
ylabel('meters or radians')
subplot(2,2,4)
plot(ts, xs(4,:), '-b','LineWidth',2)
grid on
xlabel('time(s)')
ylabel('wheel radius')
axis equal
axis xy



function [x, varargout] = uni_f(x, u, S)
%Dynamical model of the unicycle
c = cos(x(3));
s = sin(x(3));
wr = x(4);  %Wheel radius in meters

x = [x(1) + c*wr*u(1);
     x(2) + s*wr*u(1);
     x(3) + u(2);
     x(4)];

x = fix_state(x, S);

if nargout > 1
  %F-matrix
  varargout{1} = [1, 0, -s*wr*u(1), c*u(1);
                  0, 1, c*wr*u(1), s*u(1); 
                  0, 0, 1, 0;
                  0, 0, 0, 1];
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
            d(2)/r^2, -d(1)/r^2, -1,0;
            -d'/r, 0,0];
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
        H = [H;d(2)/r^2, -d(1)/r^2, -1,0];
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