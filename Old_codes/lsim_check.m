function [ ti,x] = lsim_check( )
%LSIM_CHECK Summary of this function goes here
%   Detailed explanation goes here
g2=@sense_f;
syms x
g=@(x)g2(x);

gd=diff(g(x))
gfun = matlabFunction(gd)
try
  gfun(2)
catch ME
  if (strcmp(ME.identifier ,'MATLAB:TooManyInputs'))     
   gfun()
  end
end

end


function y = sense_nl(x,S)
g=@sense_f;

y=g(x(1))*x(2);
end
function y= sense_f(x)
y=3*x+5;
end
