% Testing the fmincon for arm:

% model parameters
S.m1 = 1;
S.m2 = 1;
S.l1 = .5;
S.l2 = .5;
S.lc1 = .25;
S.lc2 = .25;
S.I1 = S.m1*S.l1^2/12;
S.I2 = S.m2*S.l2^2/12;
S.g = 9.8;

% time horizon and number of segments
tf = 2;
S.N = 50;
S.h = tf/S.N;

% Cost function parameters
S.Q = diag([1, 1, .1, .1]);
S.R = diag([.05, .05]);
S.Pf = diag([5, 5, 1, 1]);

S.f = @arm_f;
S.L = @arm_L;
S.mu = 0;

% initial state 
S.x0 = [0; 0; 0; 0];
% final desired state
S.xN = [pi/2; -pi; 0; 0];
%initial guess for fmincon states 
x0 = zeros(6*S.N-2,1); %x = [u0; x1; u1; x2;u2; x3;u3; .... xN-1;uN-1;uN] -> 6N-2 vars
%setting options:
%
options = optimoptions('fmincon','GradConstr','on','algorithm','interior-point','Display','iter','Hessian','lbfgs','GradObj','on',...
                        'MaxIter',10,'TolCon',1,'TolX',1e-5,'PlotFcns',@optimplotconstrviolation);
x = fmincon(@(x)Cost_arm(x,S), x0, [],[],[],[],[],[],@(x)constraints(x,S),options);
% initial controls 
%us = zeros(2, S.N);

%%
clear u state;
ind1 = 1:4;
ind2 = 1:2;
u = zeros(S.N,2);
for i = 0:S.N-1
u(i+1,1) = x(6*i + 1);
u(i+1,2) = x(6*i + 2);
end
state = zeros(S.N,2);
state(1,:) = S.x0(1:2)';
for i = 0:S.N-3
    state(i+2,1) = x(2*(i+1) +4*(i) + 1);
    state(i+2,2) = x(2*(i+1) +4*(i) + 2);
end
state(S.N,1:2) = S.xN(1:2)';
pt = FWD(state',S);
figure,
plot(u);
figure,
plot(pt(1,:),pt(2,:));
axis('equal');
figure,
plot(state);
 