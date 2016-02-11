function f = lecture5_3()

% draw circle
a =0:.01:2*pi; plot(cos(a),sin(a),'-','LineWidth',2)
hold on; axis equal

N = 50; % time-steps

% generate 10 initial conditions and compute the optimal trajectories
for i=1:5
  q0 = [2; 2];
  v0 = [-1.5; 0];

  vf = [0; 0];
  b = .01*i;
  
  xs = compute(q0, v0, vf, b, N);
  plot(xs(1,:), xs(2,:), '.g-')
end


function xs = compute(q0, v0, vf, b, N)

% guess for initial conditions
c2 = [1;1]; c3=[1;1]; nu1 = 1;
c = [c2; c3; nu1; 1];

c = fsolve(@(c)coeff(c, q0, v0, vf, b), c)

c2 = c(1:2); c3 = c(3:4);
tf = c(6)

h = tf/N;
xs = zeros(4,N+1);

for i=1:N+1
  t = (i-1)*h;
  [q, v, u, du] = traj(t, q0, v0, c2, c3);
  xs(:,i) = [q; v];
end


function [q, v, u, du] = traj(t, q0, v0, c2, c3)

q = c3*t^3 + c2*t^2 + v0*t + q0;
v = 3*c3*t^2 + 2*c2*t + v0;
u = 6*c3*t + 2*c2;
du = 6*c3;


function f = coeff(c, q0, v0, vf, b)
% implicit equation for the trajectory coeffcients and constraint multiplier

c2 = c(1:2);
c3 = c(3:4);
nu1 = c(5);
tf = c(6);

[q,v,u,du] = traj(tf, q0, v0, c2, c3);

f = [2*q*nu1 - du; 
     q'*q-1;
     v - vf;
     b - u'*u/2 + du'*v];
