function [ns,ks]=getnsi(wl)
% Returns two outputs-real and imaginary part of the refractive index of Si
% wl is in nm.
load ksi_wl01
load nsi_wl01
[r,~]=size(wl);
[~,loc]=ismember(wl,wl01);

ns=nsi(loc);%ns is the real part
ks=ksi(loc);%ks is the imaginary part
if r==1
    ns=ns';
    ks=ks';
end
end