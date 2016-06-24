 function [resp,fre] = bode_lin
% clear all;
 close all;
freq_logrange=unique(round(logspace(0,log10(100))));
for i=1:length(freq_logrange)
    i
[resp(i),fre(i),sys]=freq_resp_Ac(freq_logrange(i));
end
close all;
angle_resp=unwrap(angle(resp));%wrapTo2Pi(angle(resp));
s=scatter(fre,20*log10(abs(resp)),'g','filled');
alpha(s,.5)
hold on;
%figure;
p = bodeoptions;
p.FreqUnits='Hz';
%p.PhaseWrapping ='on';
%sys_H = tf(1,[1,1]);
bode(sys,p);
 hold on;
 s1=scatter(fre,angle_resp*180.0/pi,'g','filled');
 alpha(s1,.5)
 saveas(gcf,['BodeComp_','.svg']);

 end


function [response,freq,sys] = freq_resp_Ac(freqmult)
% clear all;
close all;

 
%Sampling Paramters
dt = 0.02; % 50 Hz sampling

f0 = 0.05;

num_cycles=2;
Tfinal = num_cycles/f0;

Nfinal=Tfinal/dt;
 
%Input Signal Parameters
% Base frequency of sinusoidal inputs:
%0.05;

%freqmult = 20;
freq = [freqmult]*f0; % should be 1Hz for for f0=0.05
omega=2*pi*freq;
a = 1; %Amplitude of stimulus

t=0:dt:Tfinal;
u = a * sin(omega*t)';
[~,y,sys] = Ac_SS_Fullstate(dt,Tfinal,freq);

 

 
% Drop first cycle of u & y due to transient response:
idxtrunc = (Nfinal/num_cycles)+1:Nfinal; 
utrunc = u(idxtrunc);
ytrunc = y(idxtrunc);



 
Y_fft=fft(ytrunc);
U_fft=fft(utrunc);


% To compute index of stimulus frequency, note that f0/4 correspondes to
% the resolution the FFT. This is because our total signal (over the
% duration of idxtrunc) is 80 seconds long, i.e. T0*4. So, the stimulus
% frequency is simply given by index . T0=1/f0;

 
stimidx = freqmult * (num_cycles-1) + 1; % idx = 0 corresponds to DC.

 
response=Y_fft(stimidx)/U_fft(stimidx);

 

 
end
function [response,freq] = freq_resp(freqmult)
% clear all;
close all;

 
%Sampling Paramters
dt = 0.02; % 50 Hz sampling

f0 = 0.05;

num_cycles=3;

Tfinal = num_cycles/f0;

Nfinal=Tfinal/dt;
 
%Input Signal Parameters
% Base frequency of sinusoidal inputs:
%0.05;

%freqmult = 20;4 
freq = [freqmult]*f0; % should be 1Hz for for f0=0.05
omega=2*pi*freq;
a = 1; %Amplitude of stimulus

 
t=0:dt:Tfinal;

 A=[0 1;0 -0.1];
 B=[0;1];
 C=[-47.12 0];
 D=0;
 sys_H=ss(A,B,C,D);
% make up a linear system:
%sys_H = tf(1,[1,1]);

 
% Simulate
u = a * sin(omega*t)';
y = lsim(sys_H,u,t);

 

 
% Drop first cycle of u & y due to transient response:
idxtrunc = (Nfinal/num_cycles)+1:Nfinal; 
utrunc = u(idxtrunc);
ytrunc = y(idxtrunc);



 
Y_fft=fft(ytrunc);
U_fft=fft(utrunc);


% To compute index of stimulus frequency, note that f0/4 correspondes to
% the resolution the FFT. This is because our total signal (over the
% duration of idxtrunc) is 80 seconds long, i.e. T0*4. So, the stimulus
% frequency is simply given by index . T0=1/f0;

 
stimidx = freqmult * (num_cycles-1) + 1; % idx = 0 corresponds to DC.

 
response=Y_fft(stimidx)/U_fft(stimidx);

 

 
end
