function f = lecture2_1()
% EN530.603 Applied Optimal Control
% Lecture#2, 2013/09/09
%
% Unconstrained nonlinear optimization basics (1 of 2)
% 
% Author: Marin Kobilarov, marin(at)jhu.edu


f1 = @(x1,x2)[x1; x2]'*[1, -1; -1 4]*[x1; x2];
f2 = @(x1,x2)[x1; x2]'*[-1, 1; 1 3]*[x1; x2];
f3 = @(x1,x2)(x1-x2^2)*(x1 - 2*x2^2);

figure
ezsurfc(f1, [-5 5]);
figure
ezsurfc(f2, [-5 5]);
figure
ezsurfc(f3, [-5 5]);
