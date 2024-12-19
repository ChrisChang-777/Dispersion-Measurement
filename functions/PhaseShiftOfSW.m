function [E,freq,v] = PhaseShiftOfSW(rec,dt,offset,dx,vmin,dv,vmax,fmin,fmax)
%   Summary of this function goes here.
%   [E,freq,v] = PhaseShiftOfSW(rec,dt,offset,dx,vmin,dv,vmax,fmin,fmax)
%   Detailed explanation goes here.
%   The function is for generating the dispersive energy of surface waves.
%   The function is generally adopted in Multi-channel Analysis of Surface Waves (MASW).
%
%   IN      
%           rec: the seismic record of surface waves, the size of rec is [nt*nx].
%            dt: the sampling interval in time domain (s).
%        offset: the distance (m) of the source to the first receiver.
%            dx: the spacing (m) of the receiver.
%          vmin: the minimum scanned phase velocity (m/s).
%            dv: the interval of the scanned phase velocity (m/s). 
%          vmax: the maximum scanned phase velocity (m/s).
%          fmin: the minimum frequency (Hz) of the dispersive energy.
%          fmax: the maximum frequency (Hz) of the dispersive energy.
%
%  OUT   
%             E: the normalized dispersive energy.
%          freq: the frequency (Hz) vecotor of the normalized dispersive energy.
%             v: the phase velocity (m/s) vector of the normalized dispersive energy.
%  
%   References: 
%       	Park CB, Miller RD, Xia J. Imaging dispersion curves of surface waves 
%           on multi-channel record[C]. Technical Program with Biographies SEG, 
%           68th Annual Meeting, New Orleans, L A., 1998, 1377-1380.
%
%  Author(s): Yan Yingwei
%  Copyright: 2020-2025 
%  Revision: 1.0  Date: 5/28/2020
%
%  Academy of Opto-Electronics, China Electronic Technology Group Corporation (AOE CETC)


[nt,nx] = size(rec);      % 获得地震记录的维度
nt = 4*(2^nextpow2(nt));
rec_fft = zeros(nt,nx);   % 地震记录的频谱

v=vmax:-dv:vmin;          % 相速度扫描向量
lv = length(v);           % 相速度向量的长度

Fs = 1/dt;                % 采样率
f = Fs*(0:(nt/2))/nt;     % 频率向量

% 确定频散能量图中最小频率在频率向量中的索引
ind = find(f>fmin);       
fmin_ind  = ind(1)-1;

% 确定频散能量图中最大频率在频率向量中的索引
ind = find(f>fmax);
fmax_ind = ind(1);

freq = f(fmin_ind:fmax_ind); % 频散能量图中的频率向量
lf = length(freq);

for j=1:nx
    rec_fft(:,j) = fft(rec(:,j),nt,1); % 获得地震记录的频谱
end

rec_fft_amp = abs(rec_fft);

rec_fft_n = rec_fft./rec_fft_amp;      % 对地震记录的频谱作归一化

E = zeros(lv,lf); % 频散能量矩阵


for j=1:lf
    f_ind = fmin_ind+j-1;
    for i=1:lv
        for k=1:nx
            x=offset+(k-1)*dx;
            E(i,j) = E(i,j)+exp(1i*2*pi*f(f_ind)*x/v(i))*rec_fft_n(f_ind,k);
        end
    end
    E(:,j) = abs(E(:,j)./max(abs(E(:,j))));  % 对每个频率的频散能量作归一化
end
end

