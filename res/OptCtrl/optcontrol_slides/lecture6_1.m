syms P(t) x(t) x0 tf; 

% solve P(t)
Pt = dsolve(diff(P) == P^2 - 2*P - 1, P(tf) == 0)

% solve x(t) using:
% l = P*x
% u = -2*l
% dx = x + u  and hence dx = (1-2*P)x
dsolve(diff(x) == (1 - 2*Pt)*x, x(0) == x0)

x0 = -1;
tf = 5;
t=0:.01:tf;

subplot(1,3,1)
% copy and paste solution from above
Ps = 1 - 2^(1/2)*tanh(2^(1/2)*(t - tf + (2^(1/2)*atanh(2^(1/2)/2))/2));
plot(t, Ps, 'LineWidth',3)
xlabel('t', 'FontSize', 15); ylabel('P(t)', 'FontSize', 15)

subplot(1,3,2)
plot(t, -2*Ps.*xs, 'LineWidth', 3)
xlabel('t', 'FontSize', 15); ylabel('u(t)', 'FontSize', 15)

subplot(1,3,3)
% copy and paste solution from above
xs = x0*exp(t*(2*2^(1/2) - 1) - 2*log(tanh(atanh(2^(1/2)/2) + 2^(1/2)*t - 2^(1/2)*tf) + 1))*(tanh(2^(1/2)*tf - atanh(2^(1/2)/2)) - 1)^2;
plot(t, xs, 'LineWidth', 3)
xlabel('t', 'FontSize', 15); ylabel('x(t)', 'FontSize', 15)

