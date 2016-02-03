function f = lecture5_1()

p = [3.66; -1.86];

a=[105*pi/180; 240*pi/180];

a = fsolve(@(a)angle(a, p), a)

a=[105*pi/180; 240*pi/180];

V=.3;

q = [p' a(1)];

[t,y]=ode45(@(t,q) zermelo(t,q,V), [0 18], q, []);   

% draw ship velocity
quiver(y(1:4:end,1),y(1:4:end,2), V*cos(y(1:4:end,3)), V*sin(y(1:4:end,3)))
hold on
plot(y(:,1),y(:,2),'.-g','LineWidth',2);

% draw current
%xs = rand(
%quiver(y(1:5:end,1),y(1:5:end,1), V*cos(y(1:5:end,1)), V*sin(y(1:5:end,1)))


function f = zermelo(t,q,V)
x=q(1);
y=q(2);
a=q(3);

u = -V*y;

f = [cos(a)*V + u;
     sin(a)*V;
     cos(a)^2*V];


function f = angle(a,p)
s0 = sec(a(1));
sf = sec(a(2));
t0 = tan(a(1));
tf = tan(a(2));

f = p - [(sf*(tf - t0) - t0*(sf - s0) + asinh(tf) - asinh(t0))/2; 
         s0 - sf]; 
