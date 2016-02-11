function [ Res ] = function2b( tf )
%FUNCTION2b Finds the nonlinear equation for tf
w = sqrt(3/2);
nu = -10/(4*(tanh(w*tf))^2 - (tan(w*tf))^2 -5);
Res = (nu/(2*w))*(-8*tanh(w*tf) - 2*tan(w*tf)) + 5*tf - 15;
end

