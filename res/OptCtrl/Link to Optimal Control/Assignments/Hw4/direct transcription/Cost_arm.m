function [C,gradC] = Cost_arm(x,S)
  C = 0;
  gradC = zeros(1,6*S.N-2);
  ind1 = 1:4;
  ind2 = 1:2;
  C = C + S.h/2 * x(ind2)'*S.R*x(ind2);
  gradC(ind2) = (S.R*x(ind2))';
  for i = 1:S.N-1
      xi = x(2*i +4*(i-1) + ind1);
      ui = x(6*i + ind2);
      dx = xi - S.xN;
      C  = C + S.h/2*(dx'*S.Q*dx + ui'*S.R*ui);
      gradC(2*i +4*(i-1) + ind1) = (S.h*S.Q*dx)';
      %Lxx = S.h*S.Q;
      gradC(6*i + ind2) = (S.h*S.R*ui)';
      %Luu = S.h*S.R;
  end
    C = C + S.h/2*x(6*S.N-4+ind2)'*S.R*x(6*S.N-4+ind2);
  gradC(6*S.N - 4 + ind2) = (S.R*x(6*S.N-4+ind2))';
end