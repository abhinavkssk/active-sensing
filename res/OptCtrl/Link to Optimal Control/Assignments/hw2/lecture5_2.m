function f = lecture5_2()

% draw circle
a =0:.01:2*pi; plot(cos(a),sin(a),'-','LineWidth',3)
hold on; axis equal

N = 50; % time-steps

% generate 10 initial conditions and compute the optimal trajectories
for i=1:10  
  q0 = [3*randn; 3*randn];
  v0 = [5*randn; 5*randn];
  vf = [0; 0];
  tf = 1;
  
  xs = compute(q0,v0,vf,tf, N);
  plot(xs(1,:), xs(2,:), '.g','LineWidth',2)
end


function xs = compute(q0, v0, vf, tf, N)

% guess for initial conditions
c2 = [1;1]; c3=[1;1]; nu1 = 1;
c = [c2; c3; nu1];

c = fsolve(@(c)coeff(c, q0, v0, vf, tf), c)

c2 = c(1:2); c3 = c(3:4);

h = tf/N;
xs = zeros(4,N+1);

for i=1:N+1
  t = (i-1)*h;
  [q, v] = traj(t, q0, v0, c2, c3);
  xs(:,i) = [q; v];
end


function [q, v] = traj(t, q0, v0, c2, c3)
q = c3*t^3 + c2*t^2 + v0*t + q0;
v = 3*c3*t^2 + 2*c2*t + v0;


function f = coeff(c, q0, v0, vf, tf)
% implicit equation for the trajectory coeffcients and constraint multiplier

c2 = c(1:2);
c3 = c(3:4);
nu1 = c(5);

[q,v] = traj(tf, q0, v0, c2, c3);

f = [2*q*nu1 - 6*c3; 
     q'*q-1;
     v - vf];

