function [x,A,B] = body_f( x, S)
b = S.b;
m=S.m;
A = [ 0 1;0 -b/m];

B = [0;1/m;];
end