syms t0 t1 tf x10 x20 ;
%sequence u = -1,1
x21 = -(t1-t0) + x20;
x11 = -0.5*(t1-t0)^2 + x20*(t1-t0) + x10;
x2f = (tf-t1) + x21;
x1f = 0.5*(tf-t1)^2 + x21*(tf - t1) + x11;
f1 = x2f - x1f*(t1 - tf);
disp(f1);
t1roots = solve(f1 == 0,t1);
t1alt = solve(x2f == 0, t1);
%min time to reach origin so if tf is greater than this we dont need above
%junk:
C1 = solve(x2f ==0, t1);
x1f1 = subs(x1f,t1,C1);
tfc = solve(x1f1 ==0, tf);
% if in first quadrant then we may have to use u = -1 without switching if
% tf is not enough to get into second quadrant. 
tf2 = solve(C1 == tf,tf);
syms u t;
x(t,t0,u,x10,x20) = [u*(t-t0) + x20;u*0.5*((t-t0)^2) + x20*(t-t0) + x10];
%%
clc;
conds = {0,4,0,2,3};%{t0,tf,t1,x20,x10}
%conds = {0,3,0,2,2};%{t0,tf,t1,x20,x10}
t1check = double(subs(t1roots,{t0,tf,t1,x20,x10},conds));
conds{3} = t1check(1);
disp('t1 roots');
disp(t1check);
x2fcheck = double(subs(x2f,{t0,tf,t1,x20,x10},conds));%change t1 acc
x1fcheck = double(subs(x1f,{t0,tf,t1,x20,x10},conds));%change t1 acc
disp('x1f:');
disp(x1fcheck);
disp('x2f:');
disp(x2fcheck);
dist = sqrt(x1fcheck^2 + x2fcheck^2);
disp('dist');
disp(dist);
disp('Comparing to unoptimal');
conds{3} = conds{1};
%conds = {0,2,0,-2,3};%{t0,tf,t1,x20,x10} t1 = 0
x2fcheck = double(subs(x2f,{t0,tf,t1,x20,x10},conds));%change t1 acc
x1fcheck = double(subs(x1f,{t0,tf,t1,x20,x10},conds));%change t1 acc
disp('x1f:');
disp(x1fcheck);
disp('x2f:');
disp(x2fcheck);
dist = sqrt(x1fcheck^2 + x2fcheck^2);
disp('dist');
disp(dist);
disp('Comparing to alternate2');
t1check = double(subs(t1alt,{t0,tf,t1,x20,x10},conds));
conds{3} = t1check;
disp('t1 roots');
disp(t1check);
x2fcheck = double(subs(x2f,{t0,tf,t1,x20,x10},conds));%change t1 acc
x1fcheck = double(subs(x1f,{t0,tf,t1,x20,x10},conds));%change t1 acc
disp('x1f:');
disp(x1fcheck);
disp('x2f:');
disp(x2fcheck);
dist = sqrt(x1fcheck^2 + x2fcheck^2);
disp('dist');
disp(dist);
disp('if in 1st quad then min tf needed for switching');
disp(double(subs(tf2,{t0,tf,t1,x20,x10},conds)));
disp('min tf to reach origin using -1, +1 sequence');
disp(double(subs(tfc,{t0,tf,t1,x20,x10},conds)));
    %%
%sequence u = +1,-1
x21prime = (t1-t0) + x20;
x11prime = 0.5*(t1-t0)^2 + x20*(t1-t0) + x10;
x2fprime = -(tf-t1) + x21prime;
x1fprime = -0.5*(tf-t1)^2 + x21prime*(tf - t1) + x11prime;
f1prime = x2fprime - x1fprime*(t1 - tf);
disp(f1prime);
t2roots = solve(f1prime == 0,t1);
t1altprime = solve(x2fprime == 0, t1);
t2alt = solve(x2fprime == 0, t1);
%min time to reach origin so if tf is greater than this we dont need above
%junk:
C2 = solve(x2fprime ==0, t1);
x1f2 = subs(x1fprime,t1,C2);
tfc2 = solve(x1f2 ==0, tf);
% if in third quadrant then we may have to use u = 1 without switching if
% tf is not enough to get into fourth quadrant. 
tf2 = solve(C1 == tf,tf);
%%
clc;
conds = {0,4,0,-2,-3};%{t0,tf,t1,x20,x10}
t2check = double(subs(t2roots,{t0,tf,t1,x20,x10},conds));
conds{3} = t2check(1);
disp('t1 roots');
disp(t2check);
x2fprimecheck = double(subs(x2fprime,{t0,tf,t1,x20,x10},conds));%change t1 acc
x1fprimecheck = double(subs(x1fprime,{t0,tf,t1,x20,x10},conds));%change t1 acc
disp('x1f:');
disp(x1fprimecheck);
disp('x2f:');
disp(x2fprimecheck);
dist2 = sqrt(x1fprimecheck^2 + x2fprimecheck^2);
disp('dist');
disp(dist2);
disp('Comparing to unoptimal');
conds{3} = conds{2};
%conds = {0,2,0,-2,3};%{t0,tf,t1,x20,x10} t1 = 0
x2fprimecheck = double(subs(x2fprime,{t0,tf,t1,x20,x10},conds));%change t1 acc
x1fprimecheck = double(subs(x1fprime,{t0,tf,t1,x20,x10},conds));%change t1 acc
disp('x1f:');
disp(x1fprimecheck);
disp('x2f:');
disp(x2fprimecheck);
dist2 = sqrt(x1fprimecheck^2 + x2fprimecheck^2);
disp('dist');
disp(dist2);
disp('Comparing to Alternate');
conds{3} = double(subs(t1alt,{t0,tf,t1,x20,x10},conds));
%conds = {0,2,0,-2,3};%{t0,tf,t1,x20,x10} t1 = 0
x2fprimecheck = double(subs(x2fprime,{t0,tf,t1,x20,x10},conds));%change t1 acc
x1fprimecheck = double(subs(x1fprime,{t0,tf,t1,x20,x10},conds));%change t1 acc
disp('x1f:');
disp(x1fprimecheck);
disp('x2f:');
disp(x2fprimecheck);
dist2 = sqrt(x1fprimecheck^2 + x2fprimecheck^2);
disp('dist');
disp(dist2);
disp('min tf to reach origin using 1, -1 sequence');
disp(double(subs(tfc2,{t0,tf,t1,x20,x10},conds)));
%% Plotting dynamics:
% choose t0 = 0;
%Plotting switching curves:
X = matlabFunction(x);
t = 0:0.1:3;
figure(1), clf;
X1 = X(t,0,1,4.5,-3);%(t,t0,u,x10,x20)
plot(X1(1,:),X1(2,:));
hold on;
plot(0,0,'r*');
X1 = X(t,0,-1,-4.5,3);%(t,t0,u,x10,x20)
plot(X1(1,:),X1(2,:));
grid on;
% Quadrant 2 above switching curve with t < min tf
% min tf = 2.47 for x20  = -2 and x10 = 3; 
%opt seq = -1,1
t1 = 0.1721;
tf = 1.5;
t = 0:0.1:t1;
X1 = X(t,0,-1,3,-2);%(t,t0,u,x10,x20)
plot(X1(1,:),X1(2,:),'r');
 t = t1:0.1:tf;
 X1 = X(t,t1,1,X1(2,end),X1(1,end));
 plot(X1(1,:),X1(2,:),'g');
 %Quadrant 1 and 4 above switching curve with large and small tf
 %posn(2,3) seq = {-1,1} or {-1}
 t1 = 3.6117;
 tf = 4;
 t = 0:0.1:t1;
 X1 = X(t,0,-1,3,2);%(t,t0,u,x10,x20)
 plot(X1(1,:),X1(2,:),'r');
 t = t1:0.1:tf;
 X1 = X(t,t1,1,X1(2,end),X1(1,end));
 plot(X1(1,:),X1(2,:),'g');
%similar positions in Quadrant 2 and 4 are plotted now:
%posn(-2,-3) seq = {1,-1} or {1}
 t1 = 3.6117;
 tf = 4;
 t = 0:0.1:t1;
 X1 = X(t,0,1,-3,-2);%(t,t0,u,x10,x20)
 plot(X1(1,:),X1(2,:),'r');
 t = t1:0.1:tf;
 X1 = X(t,t1,-1,X1(2,end),X1(1,end));
 plot(X1(1,:),X1(2,:),'g');
 %Quadrant 4 below switching curve:
%posn(2,-3) seq = {1,-1} or {1} (on switching curve)
 t1 = 0.1721;
 tf = 1.5;
 t = 0:0.1:t1;
 X1 = X(t,0,1,-3,2);%(t,t0,u,x10,x20)
 plot(X1(1,:),X1(2,:),'r');
 t = t1:0.1:tf;
 X1 = X(t,t1,-1,X1(2,end),X1(1,end));
 plot(X1(1,:),X1(2,:),'g');
