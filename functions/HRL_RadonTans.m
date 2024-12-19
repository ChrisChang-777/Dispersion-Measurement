% High Resolution Linear Radon Transform (HRLRT) method
%
% Usage
%   [E,f,v] = HRLRT_func(uxt,df,dt,np,vmin,vmax,fmin,fmax,x,FK)
%
% INPUT:
%   uxt, 2D seismic matrix [npts,ntrace]
%   df, Frequency sampling interval.
%   x, 1D offset info [ntrace]
%   dt, time step
%   fmin, interested frequency range minF
%   fmax, interested frequency range maxF
%   vmin, interested velocity range minV
%   vmax, interested velocity range maxV
%   FK,  A flag bit used to indicate whether to apply the F-K filter

% Auther: Chang Weishuai
% date: 2024/12/9


function [E,f,v]=HRL_RadonTans(uxt,df,dt,vmin,vmax,dv,fmin,fmax,x,FK)

if FK==1% normalization
    [uxt] = fk_filter(uxt, dt, ds);
end
np = vmax-vmin+1;
ccn=fix(1./df./dt);
D=fft(uxt,ccn);
D=D.';
lf=round(fmin./df)+1;
nf=round(fmax./df)+1;
pp=1./(vmin:vmax)';
ll0=1i.*2.*pi.*df.*(pp*x);


E_temp=zeros(np,nf);
for j=lf:nf
    l=exp(ll0.*(j-1));
    E_temp(:,j)=(l*D(:,j));
end
E=abs(E_temp);


for i=lf:nf
    E(:,i)=E(:,i)./max(E(:,i));
end
E = E(:,lf:nf);
f = fmin:df:fmax;
v = vmin:dv:vmax;

