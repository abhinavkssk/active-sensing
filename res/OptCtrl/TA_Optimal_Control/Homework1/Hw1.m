function [] = Hw1(  )
%HW1 Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%% Part A %%%%%%%%%%%%%%%
%Variables
syms x1 x2;
x = [x1;x2];

%Choose function
L = (1-x1)^2 + 100*(x2-x1^2)^2;
%Gradient and Hessian
gradL = jacobian(L,x)';
hessL = jacobian(gradL,x)';

%Matlab handles for values:
Lnum = matlabFunction(L,'vars',{x});
gradLnum = matlabFunction(gradL,'vars',{x});
hessLnum = matlabFunction(hessL,'vars',{x});

disp(['L = ' , char(L)]);
disp(['gradient L = ']); disp(gradL);
disp(['hessian L = ']); disp(hessL);



disp('Gradient descent');
%Choose step size
step = 0.002;
disp(['step = ',num2str(step)]);
%Initialpoint 
itpoints = zeros(2,1000);
x0 = [0;0];
itpoints(:,1) = x0;
disp('Starting point = '); disp(x0);
%Finding the stationary points:
i = 1;
gradcurrent = gradLnum(x0);
%figure,clf,  hold on;
while ((i < 10000) && (norm(gradcurrent) > 1e-4))
    gradcurrent = gradLnum(x0);
    x0 = x0 - step*gradcurrent; 
    disp(['iteration: ', num2str(i)]);
    %disp(['norm of gradient: ',num2str(norm(gradcurrent))]);
    %disp('gradient');
    %disp(gradcurrent);
    %plot
    %plot(i,norm(gradcurrent),'*');
    %pause(0.05);
    i = i+1;
    itpoints(:,i) = x0;
end
disp(['Final point found after',num2str(i),'is:']);
disp(x0);
disp('Value at final point:');
disp(Lnum(x0));
disp('final gradient');
disp(gradLnum(x0));
disp('hessian @ final point');
disp(hessLnum(x0));
disp(eig(hessLnum(x0)));
figure;
ezsurfc(L,[-2,2]);
hold on, plot3(itpoints(1,:),itpoints(2,:),Lnum(itpoints),'g*-');
figure;
ezcontour(L,[0,2],100);
hold on, plot(itpoints(1,:),itpoints(2,:),'mx-');
disp('Please press any key to continue');
pause;

%%%%%%%%%%%%%%%% Part B %%%%%%%%%%%%%%%
clear;
close all;
figure;

syms x u;
L = x^2 + 20*u^2;
f = x - 2*u + 3;

ezcontour(L,[-5,5],100);
hold on, ezplot(f,[-5,5]);


options = optimset([],'GradObj','on','GradConstr','on');
options.Display = 'iter';
%options.PlotFcns = @optimplotx;
y0 = [0; 2];
objfun(y0);
[yfinal] = fmincon(@objfun,y0,[],[],[],[],[],[], ...
                                   @noncon,options);
disp('yfinal'); disp(yfinal);
figure;
%plot the function:
syms x u;
y = [x;u];

L = x^2 + 20*u^2;

ezsurfc(L,[-2,2]);
hold on;
f = x - 2*u + 3;
ezplot(f,[-2,2]);
end


function [Lval,gradLval] = objfun(y1)

    syms x u;
    y = [x;u];

    L = x^2 + 20*u^2;
    gradL = jacobian(L,y)';

    Lnum = matlabFunction(L,'vars',{y});
    Lval = Lnum(y1);
    plot(y1(2),y1(1),'b*-');%plot on same figure the points 
    if(nargout >1)
        gradLnum = matlabFunction(gradL,'vars',{y});
        gradLval = gradLnum(y1);
    end
end

function [c,ceq,gc,gceq] = noncon(y1)

    syms x u;
    y = [x;u];

    f = x - 2*u + 3;
    gradf = jacobian(f,y)';

    fnum = matlabFunction(f,'vars',{y});
    c = fnum(y1);
    
    if(nargout >1)
     ceq = [];   
    end
    if(nargout > 2)
        gradfnum = matlabFunction(gradf,'vars',{y});
        gc = gradfnum(y1);
        gceq = [];
    end
end


