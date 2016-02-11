function [ ] = hw2(  )
%HW3 Summary of this function goes here
%   Detailed explanation goes here
% In the given problem:  xdot = Ax + Bu
%  A = [0, 1; 3, -1] B = [0;1]
%  Q = diag(1, 1/2) (not unique actually , R = 1/2 
% ( Q is positive semidefinite, R symm pos def
% Optimality for the cost function ensures that,
% lambdadot = -Qx - A'lambda;
% u = -Rinv*B'*lambda
% lamdba = P*x
% Pdot = -A'*P - PA + PBRinv*B'*P - Q and P(t_f) = 0
%this gives control law u

%constants for the problem: A,B,Q,R,tf
close all;
A = [0 1; 3 -1];
B = [0;1];
Q = [1 0;0 0.5];
R = 0.5;
tf = 20;
Pf = zeros(2,2);
deltat = 0.1;
%Solve riccati equation for tf = 20, Pf = 0
[t,P] = ode45(@(t,P)Pdot(t,P,A,B,Q,R),tf:-deltat:0,Pf);
m1 = figure;
for i = 1:4
    subplot(2,2,i),  plot(t,P(:,i));
    ylabel(strcat('P_',int2str(i),'(t)'));
    xlabel('time(t)');
end
%plot2svg('riccati.svg',m1);
%exportfig('pic1.pdf',...
   % 'width',3.7,...
 %   'color','rgb');
%Solve for control input and state x(t)
N = length(t);
X = zeros(2,N);
U = zeros(N,1);
X(:,N) = [-5; 5];%X(0)
for count = N:-1:2
    [xdotcount,U(count)] = xdot(t(count),X(:,count),A,B,P(count,:),R);
    X(:,count-1) = X(:,count) + deltat*xdotcount;
end
m2 = figure;
for i = 1:2
    subplot(2,2,i), plot(t,X(i,:));
    ylabel(strcat('X_',int2str(i),'(t)'));
    xlabel('time(t)');
end
subplot(2,1,2), plot(t,U);
ylabel('U(t)');
xlabel('time(t)');
%exportfig('pic2.pdf',...
 %   'width',3.7,...
%    'color','rgb');
%plot2svg('XandU.svg',m2);
end

function out1 = Pdot(t,P,A,B,Q,R)
P1 = reshape(P,2,2);
out = -(A')*P1 - P1*A + P1*B*(B')*P1*(1/R) - Q; %R is scalar 
out1 = out(:);
end

function [out1, u] = xdot(t,x,A,B,P,R)
Pt = reshape(P',2,2);
u = -(1/R)*(B')*Pt*x;
out1 = A*x + B*u;
end


