%% Surface Wave Dispersion Measurement
% This is a script for testing various dispersion measurement methods.
% Please note that the code is a parallel version,
% and it should be run on a computer that supports parallel computing.
% Auther: Chang Weishuai

% Parameter:
% fmin      - Minimum frequency of velocity spectra
% fmax      - Maximum frequency of velocity spectra
% vmin      - Minimum velocity of velocity spectra
% vmax      - Maximum velocity of velocity spectra
% dr        - Reciver spacing
% dv        - velocity resolution of phase velocity spectra
% T         - Recording time
% offset    - Source to 1st receiver distance
% x         - Receiver Distances from the source (Can be row or column vector)-Uniform spacing 
% Fs        - Sampling frequency
% time      - 1D time series [npts]

clear
addpath("functions\")
load SeisData.mat     %load seismic data

% set Parameter
fmin = 0.1;
fmax = 10;
vmin = 10;
vmax = 1000;
dv = 2;
offset = 0;
T = 6;
Fs = 1/dt; 
N_r = size(data,2);
x = (offset: dr : offset + (N_r-1)*dr);
time = (0:size(data)-1)*dt;

figure(1)
wiggle(x,time,data);
xlabel('Distance(m)')
ylabel('Time(s)')


 
figure(2);
set(gcf, 'Position', [50, 150, 1600, 800]);

%% F-k method
% time_pad  - zero pading
% space_pad - insert vertual receivers
time_pad  = 8;
space_pad = 2;
[freq,v,E_fk] = fk_fun(data, T, fmin, fmax, vmin, vmax, offset, dr, time_pad, space_pad);
subplot(2,4,1)
pcolor(freq,v,E_fk);
shading interp;
colormap turbo
set(gca,'YDir','normal','XAxisLocation','bottom');
title('F-K Method')
xlabel('Frequency (Hz)');
ylabel('Phase velocity (m/s)');
axis([fmin fmax vmin vmax]);

%% tau-p method
[freq,v,E_tp] = tp_fun (data, T, fmin, fmax, vmin, dv, vmax, offset, dr);
subplot(2,4,2)
pcolor(freq,v,E_tp);
shading interp;
colormap turbo
set(gca,'YDir','normal','XAxisLocation','bottom');
title('tau-p Method')
xlabel('Frequency (Hz)');
ylabel('Phase velocity (m/s)');
axis([fmin fmax vmin vmax]);


%% Refraction Microtremor (ReMi) method
%          x:   1D offset info [ntrace]
%   lrFlag:
%       0, positive direction and negative direction
%       1, positive direction, default 1
%       -1, negative direction
%     fkFlag:   optional output fk domain spectra, default 0
lrFlag = 1;
fkFlag = 0;
[E_ReMi,freq,v] = remi(data,x,dt,fmin,fmax,vmin,vmax,lrFlag,fkFlag);
subplot(2,4,3)
imagesc(freq,v,E_ReMi);
colormap turbo
set(gca,'YDir','normal','XAxisLocation','bottom');
title('Refraction Microtremor')
xlabel('Frequency (Hz)');
ylabel('Phase velocity (m/s)');
axis([fmin fmax vmin vmax]);

%% Space Autocorrelation (Spac) method
%          x:   1D offset info [ntrace]
%   normFlag:   frequency normalization 1/ 0 or not

normFlag = 1;
[E_spac,freq,v] = Fspac(data,x,time,normFlag,fmin,fmax,vmin,vmax);
subplot(2,4,4)
imagesc(freq,v,E_spac);
colormap turbo
set(gca,'YDir','normal','XAxisLocation','bottom');
title('Space Autocorrelation')
xlabel('Frequency (Hz)');
ylabel('Phase velocity (m/s)');
axis([fmin fmax vmin vmax]);

% %% High Resolution Linear Radon Transform (HRLRT) method (old version)
% maxitr_inner = 5; % max no. of inner iteration for LRT
% maxitr_outer = 3; % max no. of outer iteration for LRT
% lamda        = 1; % damping parameter for LRT
% pad          = 2;
% [freq,v,E_HRLRT] = HRLRT_fun(data, T, fmin, fmax, vmin, dv, vmax, offset, dr, ...
%     maxitr_inner, maxitr_outer, lamda, pad); 
% subplot(2,4,5)
% imagesc(freq,v,E_HRLRT);
% colormap turbo
% set(gca,'YDir','normal','XAxisLocation','bottom');
% title('HRLRT Method')
% xlabel('Frequency (Hz)');
% ylabel('Phase velocity (m/s)');
% axis([fmin fmax vmin vmax]);
% %%


%% High Resolution Linear Radon Transform (HRLRT) method
%   df, Frequency sampling interval.
%   FK,  A flag bit used to indicate whether to apply the F-K filter
df = 0.01;
FK = 0;

[E_HRLRT,freq,v]=HRL_RadonTans(data,df,dt,vmin,vmax,dv,fmin,fmax,x,FK);
subplot(2,4,5)
imagesc(freq,v,E_HRLRT);
colormap turbo
set(gca,'YDir','normal','XAxisLocation','bottom');
title('HRLRT Method')
xlabel('Frequency (Hz)');
ylabel('Phase velocity (m/s)');
axis([fmin fmax vmin vmax]);


% %% Frequency–Bessel Transform (F-J) method (Yan Version)
% %             r: it's a vector, the distance (m) of the receiver/station point 
% %                to the reference point.
% %     RangeFreq: the range of frequency (Hz) in the dispersive energy, such
% %                as [5 100].
% %  ParaVelocity: the parameter of the phase velocity (m/s), first is the 
% %                minimal, second is the maximal, third is the interval, such 
% %                as [50 800 1].
% %           mid: the identification of the integral method, default is
% %                'Primary', at this time, adopting the Wang's method, 'NCotes'
% %                means the Newton-Cotes integral.
% r = (0:size(data,2)-1)*dr;
% RangeFreq = [fmin fmax];
% ParaVelocity = [vmin vmax dv];
% mid = 'NCotes';
% [E_FJ,freq,v, ER,EI] = FJTransform(data,dt,r,RangeFreq,ParaVelocity, mid);
% subplot(2,4,6)
% imagesc(freq,v,E_FJ);
% colormap turbo
% set(gca,'YDir','normal','XAxisLocation','bottom');
% title('F-J Method')
% xlabel('Frequency (Hz)');
% ylabel('Phase velocity (m/s)');
% axis([fmin fmax vmin vmax]);

%% Modified Frequency–Bessel Transform (F-J) method
r = abs(nonzeros(tril(x' - x, -1)));   % Interstation spacings for different pairs of stations
Trace_norm  = 'no' ;         % 'yes' or 'no'     % Trace Normalizaton
Spec_Whitening = 'no';        % 'yes' or 'no'     % Spectral Whitening

[E_MFJ, freq, v] = Modified_FJ_fun(fmin,fmax,vmin,vmax,dv,data,Fs,N_r,r,Spec_Whitening,Trace_norm);
subplot(2,4,6)
imagesc(freq,v,E_MFJ');
colormap turbo
set(gca,'YDir','normal','XAxisLocation','bottom');
title('Modified F-J Method')
xlabel('Frequency (Hz)');
ylabel('Phase velocity (m/s)');
axis([fmin fmax vmin vmax]);

%%  Phase shift method
[E_PS,freq,v] = PhaseShiftOfSW(data,dt,offset,dr,vmin,dv,vmax,fmin,fmax);
subplot(2,4,7)
imagesc(freq,v,E_PS);
colormap turbo
set(gca,'YDir','normal','XAxisLocation','bottom');
title('Phase shift Method')
xlabel('Frequency (Hz)');
ylabel('Phase velocity (m/s)');
axis([fmin fmax vmin vmax]);


%% Phase weighted slantstacking method
%          x:   1D offset info [ntrace]
%   normFlag:   frequency normalization 1/ 0 or not
%    pwsFlag:   optional pws flag for imaging 1 or not 0
%       time:   1D time series [npts]
normFlag = 1;
pwsFlag = 1;
[E_pws,freq,v] = FPhaseshift(data,x,time,normFlag,fmin,fmax,vmin,vmax,1);
subplot(2,4,8)
imagesc(freq,v,E_pws);
colormap turbo
set(gca,'YDir','normal','XAxisLocation','bottom');
title('Phase weighted slantstacking')
xlabel('Frequency (Hz)');
ylabel('Phase velocity (m/s)');
axis([fmin fmax vmin vmax]);

delete(gcp('nocreate'));




