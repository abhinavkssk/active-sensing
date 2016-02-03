function rs = uni_br_nobs()
% EN530.603 computes nonlinear observability of unicycle with range-bearing
% from beacons
% M. Kobilarov , marin(at)jhu.edu

clear

% one beacon : ubosverbale
S.pbs = [5;
         5];    % beacon positions

%two beacons: observable
S.pbs = [5, 3;
         5, 6];    % beacon positions


syms px py th u1 u2 real

x = [px; py; th];
u = [u1; u2];

% unicycle with bearing-range
f = uni_f(x,u,S)
[z, H] = br_h(x, S)

nz = length(z);
nx = length(x);

dz = H*f;
l = [z; dz];

for i=1:3
  Dl = jacobian(l,x);

  A=simplify(subs(Dl, {px,py,th}, {1,2,pi/4}))

  % try symbolic rank and also with different values for controls
  rs = [  rank(A);
          rank(subs(A, {u1,u2}, {1,1}));
          rank(subs(A, {u1,u2}, {1,0}));
          rank(subs(A, {u1,u2}, {0,1}));
          rank(subs(A, {u1,u2}, {0,0}))]  
  
  if length(find(rs==nx))
    % rank is 3
    Dl
    break
  end
  
  % keep adding time-derivatives of z
  dz = jacobian(dz, x)*f;
  l = [l; dz];
end

if max(int32(rs))==3
  disp('observable')
else
  disp('unobservable')
end



function [f, varargout] = uni_f(x, u, S)
% dynamical model of the unicycle
c = cos(x(3));
s = sin(x(3));

f = [c*u(1);
     s*u(1);
     u(2)];

if nargout > 1
  % F-matrix
  varargout{1} = [0, 0, -s*u(1);
                  0, 0, c*u(1); 
                  0 0 1];
end



function [y, varargout] = br_h(x, S)
% bearing-range from multiple beacons

p = x(1:2);

y = [];
H = [];
for i=1:size(S.pbs, 2)
  pb = S.pbs(:, i); %i-th beacon
  d = pb - p;
  r = norm(d);
  
  th = atan2(d(2), d(1)) - x(3);
  y = [y; th; r];
  
  if nargout > 1
    % H-matrix
    H = [H;
         d(2)/r^2, -d(1)/r^2, -1;
         -d'/r, 0];
  end
end

if nargout > 1  
  varargout{1} = H;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% some more models %%%%%%%%%%%%

function [y, varargout] = mr_h(x, S)

p = x(1:2);
px = p(1);
py = p(2);

d = p;
r = norm(d);

th = atan2(d(2), d(1));
y = [th; r];

if nargout > 1
  % H-matrix
  varargout{1} = [d(2)/r^2, -d(1)/r^2, 0;
                  -d'/r, 0];
end



function [y, varargout] = b_h(x, S)

p = x(1:2);
px = p(1);
py = p(2);

d = S.p0 - p;
r = norm(d);

th = atan2(d(2), d(1)) - x(3);
y = th;

if nargout > 1
  % H-matrix
  varargout{1} = [d(2)/r^2, -d(1)/r^2, -1];
end


function [y, varargout] = r_h(x, S)

p = x(1:2);
px = p(1);
py = p(2);

d = S.p0 - p;
r = norm(d);

y = r;

if nargout > 1
  % H-matrix
  varargout{1} = [-d'/r, 0];
end
