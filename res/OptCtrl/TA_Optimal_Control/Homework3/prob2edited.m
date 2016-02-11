clc;
syms tf t0 x20 x10;
A = [0.5*((1/3)*(tf^3 - t0^3) - t0^2*(tf-t0)), 0.5*(tf^2 - t0^2) - t0*(tf-t0); ...
0.5*(tf^2 - t0^2)             , tf-t0];
C = inv(A)*[-x10 - x20*(tf-t0); -x20];
%%
conds = {0,10    ,-2,3}; %t0,tf,x20,x10
disp('conds');
disp(conds);
coeffs  = double(subs(C,{t0,tf,x20,x10},conds));
disp('coeffs');
disp(coeffs);
t = conds{1}:0.1:conds{2};
plot(t,coeffs'*[t;ones(size(t))]);