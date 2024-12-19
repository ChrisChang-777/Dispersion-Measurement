function [E,freq,v, ER,EI] = FJTransform(gather,dt,r,RangeFreq,ParaVelocity, mid, dataMode)
%   Summary of this function goes here.
%   [E,freq,v, ER, EI] = FJTransform(gather,dt,r,RangeFreq,ParaVelocity, mid)
%   Detailed explanation goes here.
%   The function is for generating the dispersive energy of surface waves.
%
%   IN      
%        gather: the raw or cross-correlation record of surface waves, the 
%                size of rec is [nt*nx].
%            dt: the sampling interval in time domain (s).
%             r: it's a vector, the distance (m) of the receiver/station point 
%                to the reference point.
%     RangeFreq: the range of frequency (Hz) in the dispersive energy, such
%                as [5 100].
%  ParaVelocity: the parameter of the phase velocity (m/s), first is the 
%                minimal, second is the maximal, third is the interval, such 
%                as [50 800 1].
%           mid: the identification of the integral method, default is
%                'Primary', at this time, adopting the Wang's method, 'NCotes'
%                means the Newton-Cotes integral.
%      dataMode: '0' reprensents the common multi-record of surface waves,in 
%                this case, the broadband source is satisfied; '1' represents
%                the data excited by the rammer.
%                
%
%  OUT   
%             E: the normalized dispersive energy of the absolute F-J spectrum.
%          freq: the frequency (Hz) vecotor of the normalized dispersive energy.
%             v: the phase velocity (m/s) vector of the normalized dispersive energy.
%            ER: the normalized dispersive energy of the real part of F-J spectrum.
%            EI: the normalized dispersive energy of the imag part of F-J spectrum.
%  
%  References: 
%     Wang, J., Wu, G., & Chen, X. (2019).Frequency-Bessel transform method
%     for effective imaging of higher-mode Rayleigh dispersion curves from
%     ambient seismic noise data. Journal of Geophysical Research: Solid Earth, 
%     124,3708¨C3723. https://doi.org/10.1029/2018JB016595.
%
%  Author(s): Yan Yingwei
%  Copyright: 2021-2025 
%  Revision: 1.0  Date: 2/19/2021
%
%  Department of Earth and Space Sciences, Southern University of Science 
%  and Technology (SUSTech).

if nargin==6
    dataMode = 0;
elseif nargin==5
    dataMode = 0;
    mid = 'Primary';
end

% Determine the base parameters of the dispersive energy image
Fs = 1/dt;             % Sampling Rate
fmin = RangeFreq(1);
fmax = RangeFreq(2);
vmin = ParaVelocity(1);
vmax = ParaVelocity(2);
dv = ParaVelocity(3);

% Do fft to the gather
[M,C] = size(gather);  % C is the number of trace
M = 1*(4*2^nextpow2(M));
gather_spec = fft(gather, M, 1)/M;

f =  Fs/2*linspace(0,1,M/2+1);

v=vmax:-dv:vmin;  

indvector1=find(f<=fmin);
ind1=length(indvector1);
indvector2=find(f>fmax);
ind2=length(f)-length(indvector2)+1;

% Do the Frequency-Bessel transform to the spectrum of the gather.
n = ind2-ind1+1;
m = length(v);
E = zeros(m,n);   % the matrix of the absolute F-J spectrum.
p = 5;            % the number of point in interval [r(j-1), r(j)]

if strcmp(mid,'Primary')
for j=ind1:ind2
    omega = 2*pi*f(j);  % current circular frequency (Hz)
    for i=1:m
        v_scan = v(i);  % current scaned velocity (m/s)
        for k=2:C
            ak = gather_spec(j, k-1)-r(k-1)*(gather_spec(j, k) - gather_spec(j, k-1))/(r(k)-r(k-1));
            bk = (gather_spec(j, k) - gather_spec(j, k-1))/(r(k)-r(k-1));
            ko = omega/v_scan; 
            
            % Adopting the method in Wang Jiannan's code.
            E(i, j-ind1+1) = E(i,j-ind1+1)+ak/ko*(r(k)*besselj(1, ko*r(k))-r(k-1)*besselj(1, ko*r(k-1)))+bk/ko*(r(k)*r(k)*besselj(1, ko*r(k))-r(k-1)*r(k-1)*besselj(1, ko*r(k-1)))+...
                            bk/(ko*ko)*(r(k)*besselj(0,ko*r(k))-r(k-1)*besselj(0,ko*r(k-1)));
            b01 = 0.0;
            b02 = 0.0;
            for ii=0:14
                b01 = b01+besselj(2*ii+1, ko*r(k-1));
                b02 = b02+besselj(2*ii+1, ko*r(k));
            end
            
            E(i, j-ind1+1) = E(i,j-ind1+1)-bk*(b02-b01+besselj(0, ko*r(k-1))*0.5+besselj(0, ko*r(k))*0.5)*(r(k)-r(k-1))/(ko*ko);
        end
    end
end
elseif strcmp(mid,'NCotes')
    for j=ind1:ind2
        omega = 2*pi*f(j);  % current circular frequency (Hz)
        for i=1:m
            v_scan = v(i);  % current scaned velocity (m/s)
            for k=2:C
                ak = gather_spec(j, k-1)-r(k-1)*(gather_spec(j, k) - gather_spec(j, k-1))/(r(k)-r(k-1));
                bk = (gather_spec(j, k) - gather_spec(j, k-1))/(r(k)-r(k-1));
                ko = omega/v_scan;
                % Adopting the Newton-Cotes integral formula. 
                rVec = linspace(r(k-1),r(k),p);
                y =  ko*besselj(0, ko*rVec);
                I = NewtonCotesQuadrature(r(k-1),r(k),y);   % Newton-Cotes integral
                E(i, j-ind1+1) = E(i,j-ind1+1)+ak/ko*(r(k)*besselj(1, ko*r(k))-r(k-1)*besselj(1, ko*r(k-1)))+bk/ko*(r(k)*r(k)*besselj(1, ko*r(k))-r(k-1)*r(k-1)*besselj(1, ko*r(k-1)))+...
                    bk/(ko*ko)*(r(k)*besselj(0,ko*r(k))-r(k-1)*besselj(0,ko*r(k-1)))-bk/(ko*ko*ko)*I;
            end
        end
    end
end

% Get the real and imag part, make an absolute operation on E
% normalized the result
ER = real(E);
EI = imag(E);
E = abs(E);
if dataMode==0
    for j=1:n
        E(:,j) = E(:,j)./max(E(:,j));
        ER(:,j) = ER(:,j)./max(abs(ER(:,j)));
        EI(:,j) = EI(:,j)./max(abs(EI(:,j)));
    end
elseif dataMode==1
    ER = ER./max(max(abs(ER)));
    EI = EI./max(max(abs(EI)));
    E = E./max(max(E));
end
E = E.^(2);
ER = abs(ER).*ER;
EI = abs(EI).*EI;
freq = f(ind1:ind2);
end

