function [c,ceq,gradc,gradceq] = constraints(x,S)
%S has system properties like mass, x0,t0, xN, N, tN
%x = [u0; x1; u1; x2;u2; x3;u3; .... xN-1;uN-1;uN] -> 6N-2 vars
%x_i = x(2*i +4*(i-1) + (1:4))
% no inequalities right now
c = [];
ceq = zeros(4*S.N,1);
ind1 = 1:4;
ind2 = 1:2;

ceq(ind1) = S.x0 - x(2+ind1) - (S.h/2)*(dynamics_arm(S.x0,x(ind2),S) + dynamics_arm(x(2+ind1),x(6 + ind2),S));
for i = 1:S.N-2
    xi = x(2*i +4*(i-1) + ind1);
    xii = x(2*(i+1) +4*(i) + ind1);%xi+1
    ui = x(6*i + ind2);
    uii = x(6*(i+1) + ind2);
    ceq(4*i + ind1) = xii - xi - (S.h/2)*(dynamics_arm(xi,ui,S) + dynamics_arm(xii,uii,S));%Ci
end
ceq(4*(S.N-1) + ind1) = S.xN - x(6*(S.N)-10 + ind1) - (S.h/2)*(dynamics_arm(x(6*(S.N)-10 + ind1),x(6*S.N - 6 + ind2),S) ...
    + dynamics_arm(S.xN,x(6*S.N-4+ind2),S));

if nargout > 2
    gradc = [];
    gradceq = zeros(4*S.N,6*S.N-2);
    gradceq(1:4,1:2) = - (S.h/2)*gradient_armu(S.x0,x(ind2),S);
    for i = 0:S.N-2 %x1->N-1 and u1->N-1
        xii = x(2*(i+1) +4*(i) + ind1);
        uii = x(6*(i+1) + ind2);
        gradceq(4*i + ind1, 6*i + 2 + ind1) = eye(4) - S.h/2 *gradient_armx(xii,uii,S);%Cixii
        gradceq(4*(i+1) + ind1, 6*i + 2 + ind1) = -eye(4) - S.h/2 *gradient_armx(xii,uii,S);%Ciixii
        gradceq(4*i + ind1, 6*(i+1) + ind2) = (-S.h/2)*gradient_armu(xii,uii,S);%Ciuii
        gradceq(4*(i+1) + ind1, 6*(i+1) + ind2) = (-S.h/2)*gradient_armu(xii,uii,S);%Ciiuii
    end
    gradceq(4*(S.N-1) + ind1,6*S.N-4 + ind2) = - (S.h/2)*gradient_armu(S.xN,x(6*S.N-4+ind2),S);
    gradceq = sparse(gradceq');
end
end
function Fx = gradient_armx(x,u,S)
    q = x(1:2);
    v = x(3:4);

    c1 = cos(q(1));
    s1 = sin(q(1));
    c2 = cos(q(2));
    s2 = sin(q(2));
    c12 = cos(q(1) + q(2));
    s12 = sin(q(1) + q(2));
    % coriolis matrix
    C = -S.m2*S.l1*S.lc2*s2*[v(2), v(1) + v(2);
                        -v(1), 0] + diag([.2;.2]);
    Cx2 = -S.m2*S.l1*S.lc2*c2*[v(2), v(1) + v(2);
                        -v(1), 0];
    Cv1 = -S.m2*S.l1*S.lc2*s2*[0, 1;
                        -1, 0];
    Cv2 = -S.m2*S.l1*S.lc2*s2*[1, 1;
                        0, 0];
    % mass elements
    m11 = S.m1*S.lc1^2 + S.m2*(S.l1^2 + S.lc2^2 + 2*S.l1*S.lc2*c2) + ...
          S.I1 + S.I2;
    m11x2 = -2*S.m2*S.l1*S.lc2*s2; %m11x1 = 0;
    
    m12 = S.m2*(S.lc2^2 + S.l1*S.lc2*c2) + S.I2;
    m12x2 = -S.m2*S.l1*S.lc2*s2;%m12x2 = 0
    
    m22 = S.m2*S.lc2^2 + S.I2;
    detM = m11*m22- m12^2;
    detMx2 = m11x2*m22 - 2*m12*m12x2;
    % mass matrix
    M = [m11, m12;
         m12, m22];  
    %Minvx1 = zeros(2,2);
    %Minv = (1/detM)*[m22, -m12; -m12, m11];
    Minvx2 = [-(detMx2*m22)/(detM^2), (-m12x2/detM + (detMx2*m12)/(detM^2));...
        (-m12x2/detM + (detMx2*m12)/(detM^2)), (m11x2/detM -(m11*detMx2)/(detM^2))];
        
     % gravity vector
    fg = [(S.m1*S.lc1 + S.m2*S.l1)*S.g*c1 + S.m2*S.lc2*S.g*c12;
          S.m2*S.lc2*S.g*c12];
    fgx1 = [-(S.m1*S.lc1 + S.m2*S.l1)*S.g*s1 - S.m2*S.lc2*S.g*s12;
          -S.m2*S.lc2*S.g*s12];
    fgx2 = [ -S.m2*S.lc2*S.g*s12;
          -S.m2*S.lc2*S.g*s12];
      
    Fx = zeros(4,4);
    Fx(3:4,1) = M\(-fgx1);
    Fx(3:4,2) = Minvx2*(u-C*v-fg) + M\(Cx2*v - fgx2);
    Fx(1:2,3:4) = eye(2);
    Fx(3:4,3) = M\(-C*[1;0] -Cv1*v);
    Fx(3:4,4) = M\(-C*[0;1] -Cv2*v);
    % acceleration
    %f(3:4) = M\(u - C*v - fg);
    %velocity
    %f(1:2) = x(3:4);
    
    %Cx(1,3:4) = 
end
function Fu = gradient_armu(x,u,S)
    q = x(1:2);
    %v = x(3:4);

    %c1 = cos(q(1));
    c2 = cos(q(2));
    %s2 = sin(q(2));
    %c12 = cos(q(1) + q(2));
    % coriolis matrix
   % C = -S.m2*S.l1*S.lc2*s2*[v(2), v(1) + v(2);
    %                    -v(1), 0] + diag([.2;.2]);

    % mass elements
    m11 = S.m1*S.lc1^2 + S.m2*(S.l1^2 + S.lc2^2 + 2*S.l1*S.lc2*c2) + ...
          S.I1 + S.I2;

    m12 = S.m2*(S.lc2^2 + S.l1*S.lc2*c2) + S.I2;

    m22 = S.m2*S.lc2^2 + S.I2;

    % mass matrix
    M = [m11, m12;
         m12, m22];  
     % gravity vector
    %fg = [(S.m1*S.lc1 + S.m2*S.l1)*S.g*c1 + S.m2*S.lc2*S.g*c12;
     %     S.m2*S.lc2*S.g*c12];
    Fu = zeros(4,2);
    Fu(1:2,1:2) = zeros(2,2);
    Fu(3:4,1) = M\([1;0]);
    Fu(3:4,2) = M\([0;1]);
    % acceleration
    %f(3:4) = M\(u - C*v - fg);
    %velocity
    %f(1:2) = x(3:4);
end
function [f] = dynamics_arm(x,u,S)
    % x is 4x1 vector containing posn and velocity and S is the system
    % props.
    f = zeros(4,1);
    q = x(1:2);
    v = x(3:4);

    c1 = cos(q(1));
    c2 = cos(q(2));
    s2 = sin(q(2));
    c12 = cos(q(1) + q(2));
    % coriolis matrix
    C = -S.m2*S.l1*S.lc2*s2*[v(2), v(1) + v(2);
                        -v(1), 0] + diag([.2;.2]);

    % mass elements
    m11 = S.m1*S.lc1^2 + S.m2*(S.l1^2 + S.lc2^2 + 2*S.l1*S.lc2*c2) + ...
          S.I1 + S.I2;

    m12 = S.m2*(S.lc2^2 + S.l1*S.lc2*c2) + S.I2;

    m22 = S.m2*S.lc2^2 + S.I2;

    % mass matrix
    M = [m11, m12;
         m12, m22];  
     % gravity vector
    fg = [(S.m1*S.lc1 + S.m2*S.l1)*S.g*c1 + S.m2*S.lc2*S.g*c12;
          S.m2*S.lc2*S.g*c12];

    % acceleration
    f(3:4) = M\(u - C*v - fg);
    %velocity
    f(1:2) = x(3:4);
    %v = v + S.h*a;
    %x = [q + S.h*v;
         %v];
end
