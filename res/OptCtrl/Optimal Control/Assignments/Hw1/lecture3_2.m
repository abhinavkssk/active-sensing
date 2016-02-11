function f = lecture3_2()
% EN530.603 Applied Optimal Control
% Lecture#3, 2013/09/11
%
% Optimization with inequality constraints (2 of 2)
% 
% Marin Kobilarov, marin(at)jhu.edu


ezplot(@(x1,x2)(x1*x2/2 + (x1+2)^2 + (x2-2)^2/2 - 2), [-6,0,-1,7])
hold on
ezcontour(@(x1, x2)(x1*exp(-x1^2-x2^2)+(x1^2+x2^2)/20), [-6,0,-1,7])
plot(-.9727,.4685,'ro');
legend('constraint','f contours','minimum');


options = optimset([],'GradObj','on','GradConstr','on');
options.Display = 'iter';

x0 = [-2; 1];
[x,fval,exitflag,output] = fmincon(@onehump,x0,[],[],[],[],[],[], ...
                                   @tiltellipse,options)


function [L, dL] = onehump(x)
% standard Matlab example

L = x(1)*exp(-x(1)^2-x(2)^2)+(x(1)^2+x(2)^2)/20;

if nargout > 1
  dL = [(1-2*x(1)^2)*exp(-x(1)^2-x(2)^2)+x(1)/10;
        -2*x(1)*x(2)*exp(-x(1)^2-x(2)^2)+x(2)/10];
end

function [c,ceq,gc,gceq] = tiltellipse(x)
% standard Matlab example

c = x(1)*x(2)/2 + (x(1)+2)^2 + (x(2)-2)^2/2 - 2;

if nargout > 1
  ceq = [];
end

if nargout > 2
  gc = [x(2)/2+2*(x(1)+2);
       x(1)/2+x(2)-2];
  gceq = [];
end