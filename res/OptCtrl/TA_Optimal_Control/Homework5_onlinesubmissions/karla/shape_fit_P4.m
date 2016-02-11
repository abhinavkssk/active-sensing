function f = shape_fit_P4
%% HW 5: To run this code you need freezeColors.m and shape_fit.m
% Farshid Alambeigi
% Problem 4
%%
close all
clear all
% static batch estimation of a shape defined as a quadratic
% function z = f(p,q) + v, and parametrized using a vector x

% workspace is the square [-s,s]x[-s,s]
s = 10;

% true shape parameter (i.e. a symmetric cup)
x_true = [1; 1; 0; 0; 0; 0];

% plot true

% measurement standard dev
std = 20;

% #of measurements
k = 8;

% generate random measurements
p = 4*s*(rand(k,1) - .5)
q = 4*s*(rand(k,1) - .5)

x = [1.2 1.3 10 1 10 1]';  %initial state
P = diag([16 16 16 16 16 16]);   %Covariance


% estimate optimal parameters x
S.R = diag(repmat(std^2, k, 1));
z = shape(p, q, x_true) + randn(k,1)*std;

xs(:, 1) = x;
zs(:, 1) = zeros(2,1);
i=1;
for n=1:2:7
  i  
 H = shape_basis(p(n:n+1), q(n:n+1))
% P=inv(H'*inv(R)*H
% x = inv(H'*inv(R)*H)*H'*inv(R)*z
% z = shape(p(n:n+1), q(n:n+1), x_true) + randn(2,1)*std
zz = z(n:n+1)
P = P - P*H'*inv(H*P*H' + S.R(1:2,1:2))*H*P;
K = P*H'*inv(S.R(1:2,1:2))
x
e = zz - H*x
x = x + K*e
xs(:,i+1) = x;
zs(:,i+1) = zz;
i=i+1;

end
%  figure(2)

% plot(xs')
% hold on
% plot(xs(1,:), xs(2,:), '-b','LineWidth',3)
% legend('true', 'estimated')

% plot estimated
figure(1)
gt = ezsurf(@(p,q)shape(p, q, x_true),[-s,s]);
colormap autumn
freezeColors

alpha(gt, 0.3)

hold on

% figure(1)
ge = ezsurf(@(p,q)shape(p,q,x),[-s,s]);
% surf(xdata1,ydata1,zdata1,'Parent',axes1,'FaceAlpha',0.3,...
%     'DisplayName','true');
colormap winter
freezeColors;

alpha(ge, .8)
shape_fit(p,q,z)
% legend('true','estimated1')
function f = shape_basis(p, q)
% quadratic function, although could be any shape
f = [p.^2, q.^2, p.*q, p, q, ones(size(p))];

function z = shape(p, q, x)
z = shape_basis(p, q)*x;

