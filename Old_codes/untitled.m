u=t;
C=sense_lin(t,x,S);
xt=x';
for i=1:size(t,1)
y_sys(i)=C(i,:)*xt(:,i);
y_mod(i)=y_sys(i)*sin(S.omega*t(i));
end
    
       


[ze,pe,ke] = ellip(5,3,30,S.omega/10,'s');
[be,ae] = zp2tf(ze,pe,ke);
H1=tf(be,ae);

tq = 0:0.001:10;
vq = interp1(t,y_mod,tq);
y_filter=lsim(H1,vq,tq);

H=-tf([3*S.omega],[2*S.m 2*S.b 0])
y_sim=lsim(H*H1,u,t);

plot(t,y_sim,'.',0:0.001:10,y_filter,'-')