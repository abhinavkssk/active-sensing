function pt = FWD(x,S)
%assuming x is a 4xn matrix
pt = [cos(x(1,:))*S.l1 + cos(x(1,:)+ x(2,:))*S.l2; sin(x(1,:))*S.l1 + sin(x(1,:)+ x(2,:))*S.l2];
end