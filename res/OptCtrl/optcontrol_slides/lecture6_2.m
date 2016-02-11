function f = lecture6_2

Pf = diag([0;0]);

A = [0 1;
     2 -1];

B = [0;
     1];

Q = diag([1; 1/2]);

R = 1/4;

tf = 5;

dt = .1;

cf = P2c(Pf);

[ts, cs] = ode45(@(tf,cf) Riccati(tf, cf, ...
                                   A, B, Q, R), [tf:-dt:0], cf)
cs = flipdim(cs,1);
ts = flipdim(ts,1);

N = size(cs,1);

subplot(1,3,1)
plot(ts, cs(:,1),ts, cs(:,2),ts, cs(:,3),  'LineWidth',3)
legend('P_{11}','P_{12}','P_{22}');
xlabel('t', 'FontSize', 15); ylabel('P(t)', 'FontSize', 15)


xs = zeros(2,N);
us = zeros(2,N);

x0 = [-4; 4];
xs(:,1) = x0;

for i=1:N,
  P = c2P(cs(i,:));
  l = P*xs(:,i);
  u = -inv(R)*B'*l;
  K = A - B*inv(R)*B'*P;
  if (i < N)
    xs(:,i+1) = xs(:,i) + dt*K*xs(:,i);
  end
  us(:,i) = u;  
end

subplot(1,3,2)
plot(ts, us, 'LineWidth', 3)
xlabel('t', 'FontSize', 15); ylabel('u(t)', 'FontSize', 15)
legend('u');


subplot(1,3,3)
plot(ts, xs(1,:), ts, xs(2,:), 'LineWidth', 3)
legend('x_1','x_2')
xlabel('t', 'FontSize', 15); ylabel('x(t)', 'FontSize', 15)


function dc = Riccati(t, c, A, B, Q, R)

P = c2P(c);
dP = -P*A - A'*P - Q + P*B*inv(R)*B'*P;
dc = P2c(dP);


function c = P2c(P)
% convert from a 2x2 symmetric matrix to a 3x1 vector of its unique entries
c = [P(1,1); P(1,2); P(2,2)];

function P = c2P(c)
% the reverse of P2c
P = [c(1), c(2);
     c(2), c(3)];