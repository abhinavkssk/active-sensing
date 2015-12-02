%function [t,x,y_filter,u,delu,S] = Ac_SS_Fullstate()
clear all;
close all;
tfinal=5;
S.k=1;

S.m=1;
S.b=0.1;
S.freq=5;
S.omega=2*pi*S.freq;
S.g=@sense_f;
S.ff=0.5;%Filter Factor



% initial state
x0 = [1;0];
S.x0 = x0;

xd= [3;0];
S.xd=xd;


options = odeset('RelTol',1e-4,'AbsTol',1e-4*ones(1,34));
options_sim = odeset('RelTol',1e-4,'AbsTol',1e-4*ones(1,4));




S.lqrN=[0;0];
S.lqrQ = diag([2; 2]);

S.lqrR = 10;
[~,S.A_sys,S.B_sys]=body_f(x0,S);
S.C_sys=[-47.1239 0];
%%
%%For Multi
S.C_sys=repmat(S.C_sys,3,1);
%%

S.D_sys=0;

[S.lqrgain,~,~] = lqr(S.A_sys,S.B_sys,S.lqrQ,S.lqrR,S.lqrN);

S.kfR = 1e-4;
S.kfQ=1e-4;
S.kfN=0;

[ze,pe,ke] = ellip(5,3,30,S.omega/10,'s');
[be,ae] = zp2tf(ze,pe,ke);
[A_f,B_f,C_f,~] = tf2ss(be,ae);

S.A_del_3filter=[S.A_sys zeros(2,5) zeros(2,5) zeros(2,5);B_f*S.C_sys(1,:) A_f A_f*0 A_f*0;...
                  B_f*S.C_sys(1,:) A_f*0 A_f  A_f*0;B_f*S.C_sys(1,:) A_f*0 A_f*0  A_f  ];
S.B_del_3filter=[S.B_sys ;0*B_f;0*B_f;0*B_f];
S.C_del_3filter=[0*S.C_sys(1,:) C_f 0*C_f 0*C_f;0*S.C_sys(1,:) 0*C_f C_f 0*C_f;0*S.C_sys(1,:) 0*C_f 0*C_f C_f];
sys=ss(S.A_del_3filter,S.B_del_3filter,S.C_del_3filter,S.D_sys);
[~,S.kfgain,~] = kalman(sys,S.kfQ,diag([S.kfR S.kfR S.kfR]),S.kfN);

sys=ss(S.A_sys,S.B_sys,S.C_sys,S.D_sys);
[~,S.kfgain_sim,~] = kalman(sys,S.kfQ,diag([S.kfR S.kfR S.kfR]),S.kfN);


[t,x]=ode45(@(t,x) Multi_SS_overall_f( t,x, S), 0:0.001:tfinal, [S.x0;S.x0-[1;0];zeros(5,1);zeros(5,1);zeros(5,1);zeros(5,1);zeros(5,1);zeros(5,1)], options); 


[t_sim,x_sim]=ode45(@(tt,xx) Multi_sim_f( tt,xx, S), 0:0.001:tfinal, [S.x0-[1;0];S.x0-[1;0]], options_sim); 

for i=1:size(t_sim,1)
delu_sim(i)=-S.lqrgain*(x_sim(i,3:4)'-S.xd);
end


for i=1:size(t,1)
delu(i)=-S.lqrgain*(x(i,3:4)'-S.xd);
end
S.ustar=-S.m*S.omega^2*cos(S.omega*t)-S.b*S.omega*sin(S.omega*t);
u=delu+S.ustar';
x=x(:,1:4);
%x=x+repmat([S.xd;S.xd]',size(t,1),1);
xstar=[cos(S.omega*t) -S.omega*sin(S.omega*t)];
x1star=cos(S.omega*t);

plot(t,x(:,3),t,x(:,1)-x1star);
legend('Est delx ','Actual delx')
title('Position x1 Vs t')
xlabel('time') % x-axis label
ylabel('X1') % y-axis label
figure;


for i=1:size(t,1)
ystar(i)=sense_nl(xstar(i,:),S);
total_y(i)=sense_nl(x(i,:),S);
end

y_sys_nl=total_y-ystar;
for i=1:size(t,1)
y_mod_nl(i)=y_sys_nl(i)*sin(S.omega*t(i));
end

[ze,pe,ke] = ellip(5,3,30,S.omega/10,'s');
[be,ae] = zp2tf(ze,pe,ke);
H1=tf(be,ae);

y_filter=lsim(H1,y_mod_nl,t);

%H=-tf([3*S.omega],[2*S.m 2*S.b 0]);
%y_sim_temp=lsim(sys,delu,t,S.x0-[1;0]);
y_sim_temp=lsim(sys,delu_sim,t_sim,S.x0-[1;0]);

y_sim=lsim(H1,y_sim_temp(:,1),t);
%figure;
plot(t,y_sim,'g.',t,y_filter,'r-')
legend('Simulated signal','Actual Signal')
title('Simulates signal Vs Actual Signal')
xlabel('time') % x-axis label
ylabel('Output Signal after filtering') % y-axis label






%%

hist_axes = gca;

h_fig = ancestor(hist_axes,'figure');

% Get x/y limits of axes using axis
limits = axis(hist_axes);




% Cache the original axes position so that axes can be repositioned to
% occupy the space used by the colorstripe if nextplot clears the histogram
% axes.
original_axes_pos = get(hist_axes,'Position');

% In GUIDE, default axes units are characters. In order for axes repositiong
% to behave properly, units need to be normalized.
hist_axes_units_old = get(hist_axes,'units');
set(hist_axes,'Units','Normalized');
% Get axis position and make room for color stripe.
pos = get(hist_axes,'pos');
stripe = 0.075;
set(hist_axes,'pos',[pos(1) pos(2)+stripe*pos(4) pos(3) (1-stripe)*pos(4)])
set(hist_axes,'Units',hist_axes_units_old);

set(hist_axes,'xticklabel','')

% Create axis for stripe
stripe_axes = axes('Parent',get(hist_axes,'Parent'),...
                'Position', [pos(1) pos(2) pos(3) stripe*pos(4)]);
				 				 
limits = axis(stripe_axes);
line(limits([1 2 2 1 1]),limits([3 3 4 4 3]),...
       'LineStyle','-',...
       'Parent',stripe_axes,...
       'Color',get(stripe_axes,'XColor'));

   % Create color stripe
    n=length(x(:,1));
    binInterval = 1/n;
    xdata = [binInterval/2 1-(binInterval/2)];
    %limits(1:2) = range;

        C = (1:n)/n;

        x_range=x(:,1);
    y_range=1.5*x_range.^2+5*x_range;
        C1 = mat2gray(y_range)';
        
    % image(X,Y,C) where C is the RGB color you specify. 
    image(xdata,[0 1],repmat(C1, [1 1 3]),'Parent',stripe_axes);


%set(stripe_axes,'yticklabel','')
%axis(stripe_axes,limits);

% Put a border around the stripe.
%{
figure;
plot(t,C,t,C1)
figure;
plot(t,x_range);figure;plot(t,y_range)
%}





