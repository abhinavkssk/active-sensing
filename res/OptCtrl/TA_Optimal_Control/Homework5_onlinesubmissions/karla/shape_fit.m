function f = shape_fit(p,q,z)
%% HW5-problem 4
% Farshid Alambeigi

% close all
% static batch estimation of a shape defined as a quadratic
% function z = f(p,q) + v, and parametrized using a vector x

% workspace is the square [-s,s]x[-s,s]
s = 10;

% true shape parameter (i.e. a symmetric cup)
x_true = [1; 1; 0; 0; 0; 0];

% plot true
% figure(1)
% hold on
% gt = ezsurf(@(p,q)shape(p, q, x_true),[-s,s]);
% colormap pink
% freezeColors
% alpha(gt, 0.3)

% measurement standard dev
std = 20;

% #of measurements
k = 8;

% generate random measurements
% p = 4*s*(rand(k,1) - .5)
% q = 4*s*(rand(k,1) - .5)
% z = shape(p, q, x_true) + randn(k,1)*std

% estimate optimal parameters x
R = diag(repmat(std^2, k, 1));
H = shape_basis(p, q);
x = inv(H'*inv(R)*H)*H'*inv(R)*z

% plot estimated
hold on
% figure(1)
ge = ezsurf(@(p,q)shape(p,q,x),[-s,s]);
colormap gray
freezeColors
alpha(ge, .8)
legend('true', 'estimated','batch estimation')


function f = shape_basis(p, q)
% quadratic function, although could be any shape
f = [p.^2, q.^2, p.*q, p, q, ones(size(p))];

function z = shape(p, q, x)
z = shape_basis(p, q)*x;

