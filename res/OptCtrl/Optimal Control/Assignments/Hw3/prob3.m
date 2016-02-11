 p0 = sym('p0',[2,1]);
 p0 = sym(p0,'real');
 t = sym('t',[4,1]);
 t = sym(t,'real');
 l = sym('l',[2,1]);
 l = sym(l,'real');
 %v = sym('v',[2,1]);
 %v = sym(v,'real');
 u = sym('u',[2,1]);
 u = sym(u,'real');
 p = [cos(t(1))*l(1) + cos(t(1) + t(2))*l(2);...
     cos(t(1))*l(1) + cos(t(1) + t(2))*l(2)];
 r0 = sym('r0','real');
 C(t) = r0^2 - (p-p0)'*(p-p0);
 %%
 gradC(t) = gradient(C,t);
 %disp(gradC);
 cdot = gradC' *[t(3);t(4);u(1);u(2)];
 disp(cdot);
 gradCdot(t) = gradient(cdot,t);
 cddot = gradCdot' * [t(3);t(4);u(1);u(2)];
 disp(cddot);