
clc;
clear all;
close all;
n =5;
f = 5;

[ze,pe,ke] = ellip(n,3,30,2*pi*f/10,'s');
[be,ae] = zp2tf(ze,pe,ke);
H=tf(be,ae)
t=0:1e-4:10;

u=10*t+sin(2*pi*f*t);
plot(t,u,'g');
figure
y=lsim(H,u,t);
plot(t,y,'b');
