function [xs,us] = ddp_armmod(varargin)
%figures clear;
figure(1), clf;
figure(2), clf;
% model parameters
S.m1 = 1;
S.m2 = 1;
S.l1 = .5;
S.l2 = .5;
S.lc1 = .25;
S.lc2 = .25;
S.I1 = S.m1*(S.l1^2)/12;
S.I2 = S.m2*(S.l2^2)/12;
S.g = 9.8;

%Obstacle parameters
%  S.ck = 1;
%  S.p0 = [0.8;0.6];
%  S.rad = 0.4;
 S.mu  = 0.3;
 S.p0 = [1.6;0.3];
 S.rad = 0.0;
% time horizon and number of segments
tf = 2;
S.N = 128;
S.h = tf/S.N;

% Cost function parameters
S.Q = diag([1, 1, .1, .1]);
S.R = diag([.05, .05]);
S.Pf = diag([5, 5, 1, 1]);

S.f = @arm_f;
S.L = @arm_L;
S.mu = 1e-2;
% initial state 
x0 = [0; 0; 0; 0];
% final desired state
%S.xf = [pi/4; pi/4; 0; 0];
S.xf = [pi/2; -pi; 0; 0];
% initial controls 
%us = zeros(2, S.N);
us = randn(2,S.N);
%us(1,1:end) = -1;
%us = [repmat([.1;0.1], 1, S.N/2), repmat(-[.1;0.1], 1, S.N/2)]/5;
%us = [repmat([.1;0], 1, N/2), repmat(-[.1;0], 1, N/2)]/5;
if nargin == 1
    us = varargin{1};
end
xs = ddp_traj(x0, us, S);

Jmin = ddp_cost(xs, us,  S)
figure(1),
subplot(1,2,1)

plot(xs(1,:), xs(2,:), '-b')
hold on

plot(S.xf(1),S.xf(2),'*g')

S.a = 1;
figure(2), hold on;
%plot(S.p0,'*r');
plot(S.p0(1) + S.rad*cos(0:0.1:2*pi),S.p0(2) + S.rad*sin(0:0.1:2*pi),'r');
traj = FWD(xs,S);
plot(traj(1,:),traj(2,:),'b');
J = ddp_cost(xs, us, S)
S.ck = 5;
for j = 1:25  
    for i=1:10
      [dus, V, Vn, dV, a] = ddp(x0, us, S);

      % update controls
      us = us + dus;

      S.a = a;   %reuse old step-size for efficiency

      % update trajectory using new controls
      xs = ddp_traj(x0, us, S);
      figure(1),
      plot(xs(1,:), xs(2,:), '-b');
      figure(2),
      traj = FWD(xs,S);
      plot(traj(1,:),traj(2,:),'b');
      J = ddp_cost(xs, us, S)
    %  keyboard
    end
   % us = zeros(2, S.N);
%      if j < 5
%         S.ck = 3*S.ck;
%      else
%          S.ck =  1.4*S.ck;
%      end
     S.ck = (1 + 5/j)*S.ck;
     %S.mu = 10*S.mu;%increase regularization
     S.a = 0.2;%reinitialize step size
     figure(1), clf, hold on;
     axis('equal');
   figure(2), clf, hold on;
   axis('equal');
   plot(S.p0(1) + S.rad*cos(0:0.1:2*pi),S.p0(2) + S.rad*sin(0:0.1:2*pi),'r');
    traj = FWD(xs,S);
    plot(traj(1,:),traj(2,:),'g');
 end
figure(1),
subplot(1,2,1),
plot(xs(1,:), xs(2,:), '-m');

J = ddp_cost(xs, us, S)

xlabel('q_1')
ylabel('q_2')

subplot(1,2,2)

plot(0:S.h:tf-S.h, us(1,:),0:S.h:tf-S.h, us(2,:));
xlabel('sec.')
ylabel('Nm')
legend('u_1','u_2')

 
function pt = FWD(x,S)
%assuming x is a 4xn matrix
pt = [cos(x(1,:))*S.l1 + cos(x(1,:)+ x(2,:))*S.l2; sin(x(1,:))*S.l1 + sin(x(1,:)+ x(2,:))*S.l2];


function [L, Lx, Lxx, Lu, Luu] = arm_L(k, x, u, S)
% arm cost function with penalized term

dx = x - S.xf;

dx0 = FWD(x,S) - S.p0;%distance vec from obstacle to current position
penality = dx0'*dx0-S.rad^2;
if (k==S.N+1)
  L = dx'*S.Pf*dx/2;
  Lx = S.Pf*dx;
  Lxx = S.Pf;
  Lu = [];
  Luu = [];
else
  L = S.h/2*(dx'*S.Q*dx + u'*S.R*u);
  if penality < 0
      disp('obstacle observed');
      Ptx = [-sin(x(1))*S.l1 - sin(x(1) + x(2))*S.l2,   -sin(x(1) + x(2))*S.l2, 0, 0;
             cos(x(1))*S.l1 + cos(x(1) + x(2))*S.l2,    cos(x(1) + x(2))*S.l2,  0, 0];
         
      Ptxx1 = [-cos(x(1))*S.l1 - cos(x(1)+x(2))*S.l2,   - cos(x(1)+x(2))*S.l2,  0,    0;
                -cos(x(1)+x(2))*S.l2,                   -cos(x(1)+x(2))*S.l2,   0,    0;
                0,                                      0,                      0,    0;
                0,                                      0,                      0,    0];
            
     Ptxx2 = [-sin(x(1))*S.l1 - sin(x(1)+x(2))*S.l2,    - sin(x(1)+x(2))*S.l2,  0,    0;
                -sin(x(1)+x(2))*S.l2,                   -sin(x(1)+x(2))*S.l2,   0,    0;
                0,                                      0,                      0,    0;
                0,                                      0,                      0,    0];
      L = L + S.h*S.ck*penality^2;      
      Lx = S.h*S.Q*dx + 4*S.h*S.ck*penality*Ptx'*dx0;
      Lxx = S.h*S.Q + 4*S.h*S.ck*penality*(Ptxx1*dx0(1) + Ptxx2*dx0(2) +(Ptx')*Ptx);
      Lu = S.h*S.R*u;
      Luu = S.h*S.R;
  else
      Lx = S.h*S.Q*dx;
      Lxx = S.h*S.Q;
      Lu = S.h*S.R*u;
      Luu = S.h*S.R;
  end
end



function [x, A, B] = arm_f(k, x, u, S)
% arm discrete dynamics
% set jacobians A, B to [] if unavailable

q = x(1:2);
v = x(3:4);

c1 = cos(q(1));
c2 = cos(q(2));
s2 = sin(q(2));
c12 = cos(q(1) + q(2));

% coriolis matrix
C = -S.m2*S.l1*S.lc2*s2*[v(2), v(1) + v(2);
                    -v(1), 0] + diag([.2;.2]);

% mass elements
m11 = S.m1*S.lc1^2 + S.m2*(S.l1^2 + S.lc2^2 + 2*S.l1*S.lc2*c2) + ...
      S.I1 + S.I2;

m12 = S.m2*(S.lc2^2 + S.l1*S.lc2*c2) + S.I2;

m22 = S.m2*S.lc2^2 + S.I2;

% mass matrix
M = [m11, m12;
     m12, m22];

% gravity vector
fg = [(S.m1*S.lc1 + S.m2*S.l1)*S.g*c1 + S.m2*S.lc2*S.g*c12;
      S.m2*S.lc2*S.g*c12];


% acceleration
a = pinv(M)*(u - C*v - fg);
if norm(a) > 1e5
    disp('Very High Acceleration');
end
v = v + S.h*a;

x = [q + S.h*v;
     v];

% leave empty to use finite difference approximation
A= [];
B= [];
