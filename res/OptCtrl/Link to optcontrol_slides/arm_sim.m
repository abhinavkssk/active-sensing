function f = arm_sim()
% Example simulation of two-link arm dynamics using discrete
% dynamics
% 
% Marin Kobilarov marin(at)jhu.edu

% model parameters
S.m1 = 1;
S.m2 = 1;
S.l1 = .5;
S.l2 = .5;
S.lc1 = .25;
S.lc2 = .25;
S.I1 = S.m1*S.l1/12;
S.I2 = S.m2*S.l2/12;
S.g = 9.8;

% time horizon and number of segments
tf = 2;
N = 128;

S.h = tf/N;  %time-step

% initial state 
x0 = [0; 0; 0; 0];

% controls 
us = zeros(2, N);

% states
xs = zeros(4, N+1);
xs(:,1) = x0;

for k=1:N,
  xs(:, k+1) = arm_f(k, xs(:,k), us(:,k), S);  
end

plot(xs(1,:), xs(2,:), '-b')


function [x, A, B] = arm_f(k, x, u, S)
% arm discrete dynamics
% set jacobians A, B to [] if unavailable
%
% the following parameters should be set:
% S.m1  : mass of first body
% S.m2  : mass of second body
% S.l1  : length of first body
% S.l2  : length of second body
% S.lc1 : distance to COM
% S.lc2 : distance to COM
% S.I1  : inertia 1
% S.I2  : inertia 2
% S.g   : gravity
%
% S.h : time-step

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
a = inv(M)*(u - C*v - fg);
v = v + S.h*a;

x = [q + S.h*v;
     v];

% leave empty to use finite difference approximation
A= [];
B= [];
