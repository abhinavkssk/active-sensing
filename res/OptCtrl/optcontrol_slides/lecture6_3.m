function f = lecture6_3

Pf = diag([2;0]);

A = [0 1;
     2 -1];

B = [0;
     1];

Q = diag([2; 0]);

R = 0.005;

tf = 5;
dt = .1;

xd = [1;0];

sf = - Pf*xd;

csf = [P2c(Pf); sf];

[ts, css] = ode45(@(tf, csf) Riccati(tf, csf, ...
                                     A, B, Q, R, xd), [tf:-dt:0], csf);

css = flipdim(css,1);
ts = flipdim(ts,1);

N = size(css,1);

cs = css(:, 1:3);
ss = css(:, 4:5)';

x0 = [-4; 4];

subplot(1,3,1);
% copy and paste solution from above
plot(ts, css(:,1), ts, css(:,2), ts, css(:,3), 'LineWidth',3);
legend('P_{11}','P_{12}','P_{22}');
xlabel('t', 'FontSize', 15); ylabel('P(t)', 'FontSize', 15);

xs = zeros(2,N);
us = zeros(2,N);
xs(:,1) = x0;

for i=1:N,
  P = c2P(cs(i, :));
  s = ss(:,i);
  l = P*xs(:,i) + s;
  u = -inv(R)*B'*l;
  if (i < N)
    xs(:,i+1) = xs(:,i) + dt*(A*xs(:,i) + B*u);
  end
  us(:,i) = u;  
end

subplot(1,3,2)
plot(ts, xs(1,:), ts, xs(2,:), 'LineWidth', 3)
legend('x_1','x_2')
xlabel('t', 'FontSize', 15); ylabel('x(t)', 'FontSize', 15)

subplot(1,3,3)
plot(ts, us, 'LineWidth', 3)
xlabel('t', 'FontSize', 15); ylabel('u(t)', 'FontSize', 15)
legend('u');


function dcs = Riccati(t, cs, A, B, Q, R, xd)

P = c2P(cs(1:3));
s = cs(4:5);
dP = -P*A - A'*P - Q + P*B*inv(R)*B'*P;
ds = - (A' - P*B*inv(R)*B')*s + Q*xd;
dcs = [P2c(dP); ds];

function c = P2c(P)
c = [P(1,1); P(1,2); P(2,2)];

function P = c2P(c)
P = [c(1), c(2);
     c(2), c(3)];