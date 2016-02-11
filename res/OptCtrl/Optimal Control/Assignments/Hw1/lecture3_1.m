function f = lecture3_1()
% EN530.603 Applied Optimal Control
% Lecture#3, 2013/09/11
%
% Unconstrained nonlinear optimization basics (1 of 2)
% 
% Marin Kobilarov, marin(at)jhu.edu

q = 1/4;
r = 1;
m = -.5;
c = 1;

% cost function and its gradient
L = @(y1,y2)(q*y1^2 + r*y2^2)/2;
dL = @(y1,y2)[q*y1; r*y2];

% constraint and its gradient
f = @(y1,y2)(y1 + m*y2 - c);
df = @(y1,y2)([1; m]);

f2 = @(y2)(- m*y2 + c);

ezcontour(L, [-2 2],100);
hold on
ezplot(f, [-2 2]);

% minimum
y1s = r*c/(r + m^2*q);
y2s = m*q/r*y1s;
plot(y1s,y2s,'*g','LineWidth',2)

legend('cost contours', 'constraint', 'minimum','gradients');

% the gradient at optimum
gs = dL(y1s, y2s);
quiver(y1s, y2s, gs(1), gs(2));

% plot gradients along constraint
for y2=-2:.2:2,
  y1 = f2(y2);
  g = dL(y1, y2);
  quiver(y1, y2, g(1), g(2));
  plot(y1,y2,'o')
end

