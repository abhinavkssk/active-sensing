% Defining quantities:
dt=10^(-6); %fraction of a sec.
sn=1.5*(10^(-5));
su=3*(10^(-9));
sv=3*(10^(-6));
angrate=0.02; % rad/sec.
theta0=0;
beta0=1.7*(10^(-7));
P0=diag([10^(-4),10^(-12)]);
ci=0.95; % Confidence interval area.
T=10^(3); % Total time-steps.

F=zeros(2,2);
F(1,2)=-1;
Phi=eye(2)+dt*F;
G=[1,0]';
Gamma=dt*G;

Q=[sv^2 0; 0 su^2]*dt;
H=eye(2);
R=zeros(2,2);
R(1,1)=sn^2;

% True angle:
thetak=dt*angrate*ones(1,T+1);
thetak(1)=theta0;
thetak=cumsum(thetak,2);

% Observed angle:
zk=thetak+(sn*randn(1,T+1));

% Real beta values:
betak=randn(1,T+1)*su*dt;
betak(1)=beta0;
betak=cumsum(betak,2);

% Control values:
u=angrate*ones(1,T+1)+betak+sv*randn(1,T+1);

% Real x:
xreal=[thetak;betak];

% Observed x:
xobs=[zk;betak];

% Just for pre-allocation:
xhat=xobs; 
xpred=xhat; 
cp=zeros(1,T+1);
cm=cp;

Phat=P0;
cp(1)=1.96*sqrt(Phat(1,1));
cm(1)=-1.96*sqrt(Phat(1,1));

for i=1:T
    %Prediction step;
    xpred(:,i+1)=(Phi*xhat(:,i))+(Gamma*u(i)); 
    Ppred=(Phi*Phat*(Phi'))+Q;
    %Estimation step:
    K=(Ppred*H')/((H*Ppred*H')+R); 
    Phat=(eye(2)-(K*H))*Ppred;
    xhat(:,i+1)=xpred(:,i+1)+(K*(zk(:,i+1)-(H*xpred(:,i+1))));
    % Confidence interval storage:
    cp(i+1)=1.96*sqrt(Phat(1,1));
    cm(i+1)=-1.96*sqrt(Phat(1,1));
end

% Plotting:

red=[242/255 80/255 80/255];
green=[67/255 250/255 131/255];
yellow=[250/255 174/255 67/255];

subplot(1,2,1)
plot(0:T,xobs(1,:),'b','LineWidth',.5)
hold on
plot(0:T,xreal(1,:),'--','Color',green,'LineWidth',2)
plot(0:T,xhat(1,:),'Color',red,'LineWidth',3)
ylabel('Radians'); xlabel('Time-Steps');
legend('Observed Angle','True Angle','Estimated Angle');

subplot(1,2,2)
plot(0:T,xreal(1,:),'--','Color',green,'LineWidth',1.1)
hold on
plot(0:T,xhat(1,:),'Color',red,'LineWidth',2)
plot(1:T,xhat(1,2:end)+cp(1,2:end),'--','Color',yellow,'LineWidth',1)
plot(1:T,xhat(1,2:end)+cm(1,2:end),'--','Color',yellow,'LineWidth',1)
ylabel('Radians'); xlabel('Time-Steps');
legend('True Angle','Estimated Angle','95% Confidence Interval');

fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',13) 