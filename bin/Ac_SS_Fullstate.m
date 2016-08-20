function [t,y_filter,sys_with_filter] = Ac_SS_Fullstate(dt_sim,tfinal_bode,freq_bode)
%clear all;
close all;
addpath('../include','-end')
tfinal=40;
S.k=1;

S.m=1;
S.b= 1.767*S.m;
S.freq=2;
S.cutoffFreq=0.5;
S.cutoffFactor=S.cutoffFreq/S.freq;
S.omega=2*pi*S.freq;
S.g=@sense_f;
S.ff=0.5;%Filter Factor

dt=0.02;
S.ampl=1;
% initial state
S.x0star = [S.ampl*1;0];
S.x0 = [-3.5;0];

xd= [0.36;0];
S.xd=xd;
if nargin == 3
    S.x0 = [1;0];

    xd= [3;0];

end
S.d1=3;

options = odeset('RelTol',1e-4,'AbsTol',1e-4*ones(1,14));
options_sim = odeset('RelTol',1e-4,'AbsTol',1e-4*ones(1,4));




S.lqrN=[0;0];
S.lqrQ = diag([2; 2]);

S.lqrR = 1;
[~,S.A_sys,S.B_sys]=body_f(S.x0,S);
S.C_sys=[-S.d1*S.omega*S.ampl/(2*S.m) 0];
S.D_sys=0;

[S.lqrgain,~,~] = lqr(S.A_sys,S.B_sys,S.lqrQ,S.lqrR,S.lqrN);

S.kfR = 1e-4;
S.kfQ=1e-4;
S.kfN=0;


% [ze,pe,ke] = ellip(5,3,30,S.cutoffFactor*S.omega,'s');
% [be,ae] = zp2tf(ze,pe,ke);

[ze,pe,ke] =butter(5,S.cutoffFactor*S.omega,'s'); %ellip(5,3,30,S.cutoffFactor*S.omega,'s');
[be,ae] = zp2tf(ze,pe,ke);

S.be=be;S.ae=ae;
H1=tf(be,ae);

[A_f,B_f,C_f,~] = tf2ss(be,ae);

S.A_del_filter=[S.A_sys zeros(2,5);B_f*S.C_sys A_f];
S.B_del_filter=[S.B_sys ;0*B_f];
S.C_del_filter=[0*S.C_sys C_f];
sys=ss(S.A_del_filter,S.B_del_filter,S.C_del_filter,S.D_sys);
[~,S.kfgain,~] = kalman(sys,S.kfQ,S.kfR,S.kfN);

sys=ss(S.A_sys,S.B_sys,S.C_sys,S.D_sys);
[~,S.kfgain_sim,~] = kalman(sys,S.kfQ,S.kfR,S.kfN);

if nargin<3
    S.bode=false;
    sys_with_filter=sys*H1;
end

if nargin == 3
   S.bode=true;
   dt=dt_sim;
   tfinal=tfinal_bode;
   S.freq_bode=freq_bode;
   S.omega_bode=2*pi*S.freq_bode;
   sys_with_filter=sys*H1;
end

[t,x]=ode45(@(t,x) SS_overall_f( t,x, S), 0:dt:tfinal, [S.x0;S.x0-S.x0star;zeros(5,1);zeros(5,1)], options); 

%[t,x]=ode45(@(t,x) SS_overall_f_bode( t,x, S), 0:dt_sim:tfinal, [S.x0;S.x0-[1;0];zeros(5,1);zeros(5,1)], options); 


[t_sim,x_sim]=ode45(@(tt,xx) sim_f( tt,xx, S), 0:dt:tfinal, [S.x0-S.x0star;S.x0-S.x0star], options_sim); 

for i=1:size(t_sim,1)
delu_sim(i)=-S.lqrgain*(x_sim(i,3:4)'-S.xd);

if(S.bode)delu_sim(i)=sin(S.omega_bode*t_sim(i));end

end


for i=1:size(t,1)
delu(i)=-S.lqrgain*(x(i,3:4)'-S.xd);

if(S.bode)delu(i)=sin(S.omega_bode*t(i));end

end
S.ustar=S.ampl*(-S.m*S.omega^2*cos(S.omega*t)-S.b*S.omega*sin(S.omega*t));
u=delu+S.ustar';
x=x(:,1:4);
%x=x+repmat([S.xd;S.xd]',size(t,1),1);
xstar=S.ampl*[cos(S.omega*t) -S.omega*sin(S.omega*t)];
x1star=S.ampl*cos(S.omega*t);
x2star=-S.omega*S.ampl*sin(S.omega*t);

plot(t,x(:,3),t,x(:,1)-x1star);
legend('Est delx ','Actual delx')
legend boxoff;
title('Position vs Time')
xlabel('time (s)') % x-axis label
ylabel('\deltax_1 (m)') % y-axis label
if(S.bode)saveas(gcf,['figs/estActX_',num2str(freq_bode),'.svg']);end
figure;

plot(t,x(:,4),t,x(:,2)-x2star);
legend('Est delv ','Actual delv')
legend boxoff;
title('Velocity vs Time')
xlabel('time (s)') % x-axis label
ylabel('\deltax_2 (m/s)') % y-axis label
if(S.bode)saveas(gcf,['figs/estActvel_',num2str(freq_bode),'.svg']);end
figure;

for i=1:size(t,1)
ystar(i)=sense_nl(xstar(i,:),S);
total_y(i)=sense_nl(x(i,:),S);
end

y_sys_nl=total_y-ystar;
for i=1:size(t,1)
y_mod_nl(i)=y_sys_nl(i)*sin(S.omega*t(i));
end



y_filter=lsim(H1,y_mod_nl,t);

%H=-tf([3*S.omega],[2*S.m 2*S.b 0]);
%y_sim_temp=lsim(sys,delu,t,S.x0-S.x0star);
y_sim_temp=lsim(sys,delu_sim,t_sim,S.x0-S.x0star);

y_sim=lsim(H1,y_sim_temp,t);

%figure;
plot(t,y_sim,'g.',t,y_filter,'r-')
[lengh,objh,outh,outm]=legend('Simulated signal','Actual Signal')
%set(objh,'linewidth',10);

legend boxoff;
title('Simulated Signal vs Actual Signal')
xlabel('time (s)') % x-axis label
ylabel('Output Signal after filtering') % y-axis label
if(S.bode)saveas(gcf,['figs/Outpout_',num2str(freq_bode),'.svg']);end

if(S.bode)
    return 
end
 %set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);


% set(gca,'Ytick',[]);



%set(stripe_axes,'yticklabel','')
%axis(stripe_axes,limits);

% Put a border around the stripe.
%close all;
%figure;
%plot(t,C,t,C1)
x_range=x(:,1);

y_range=1.5*x_range.^2+5*x_range;
 C1 = mat2gray(y_range)';
figure;
%subplot(2,1,2);
plot(t,x_range)
%legend('t','x_range')
title('Position vs Time')
xlabel('time (s)') % x-axis label
ylabel('Position') % y-axis label

%  subplot(2,1,1);
%  plot(t,y_range)
%  %legend('t','x_range')
%  title('Scene Vs time')
%  xlabel('time (s)') % x-axis label
%  ylabel('Scene') % y-axis label
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%handle=subplot(212);

hist_axes =gca;

h_fig = ancestor(hist_axes,'figure');

% Get x/y limits of axes using axis
limits = axis(hist_axes);




% Cache the original axes position so that axes can be repositioned to
% occupy the space used by the colorstripe if nextplot clears the histogram
% axes.
%original_axes_pos = get(hist_axes,'Position');

% In GUIDE, default axes units are characters. In order for axes repositiong
% to behave properly, units need to be normalized.
hist_axes_units_old = get(hist_axes,'units');
set(hist_axes,'Units','Normalized');
% Get axis position and make room for color stripe.
pos = get(hist_axes,'pos');
stripe = 0.075;
%set(hist_axes,'pos',[pos(1)+stripe*pos(3) pos(2) (1-stripe)*pos(3) pos(4)])
set(hist_axes,'pos',[pos(1)+stripe*pos(3)+0.05 pos(2)+stripe*pos(4)+0.05 (1-stripe)*pos(3) (1-stripe)*pos(4)])
set(hist_axes,'Units',hist_axes_units_old);

%set(hist_axes,'XTickLabel','')

% Create axis for stripe
% stripe_axes = axes('Parent',get(hist_axes,'Parent'),...
%                 'Position', [pos(1) pos(2) stripe*pos(3) pos(4)]);
            
stripe_axes = axes('Parent',get(hist_axes,'Parent'),...
                'Position', [pos(1)+stripe*pos(3)+0.05 pos(2)-0.08 (1-stripe)*pos(3) stripe*pos(4)]);


stripe_axes1 = axes('Parent',get(hist_axes,'Parent'),...
                'Position', [pos(1)-0.08 pos(2)+stripe*pos(4)+0.05 stripe*pos(3) (1-stripe)*pos(4)]);
				 				 
limits = axis(stripe_axes);
line(limits([1 2 2 1 1]),limits([3 3 4 4 3]),...
       'LineStyle','-',...
       'Parent',stripe_axes,...
       'Color',get(stripe_axes,'XColor'));

limits1 = axis(stripe_axes1);
line(limits1([1 2 2 1 1]),limits1([3 3 4 4 3]),...
       'LineStyle','-',...
       'Parent',stripe_axes1,...
       'Color',get(stripe_axes1,'XColor'));


   % Create color stripe
    n=length(x(:,1));
    binInterval = 1/n;
    xdata = [binInterval/2 1-(binInterval/2)];
    %limits(1:2) = range;

        C = (1:n)/n;

        x_range=x(:,1);
    y_range=1.5*x_range.^2+5*x_range;
        C1 = mat2gray(y_range);
        x_range2=min(x_range):0.01:max(x_range);
         y_range2=1.5*x_range2.^2+5*x_range2;
        C2 = mat2gray(y_range2);
   
    % image(X,Y,C) where C is the RGB color you specify. 
    image(xdata,[0 1],repmat(C1', [1 1 3]),'Parent',stripe_axes);


    image(xdata,[1 0],repmat(C2', [1 1 3]),'Parent',stripe_axes1);

set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca,'Xtick',[]);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
set(gca,'Ytick',[]);

figure;
p = bodeoptions;
p.FreqUnits='Hz';
p.PhaseWrapping ='on';
%sys_H = tf(1,[1,1]);
bode(sys_with_filter,p);
%hist_axes1=hist_axes;



