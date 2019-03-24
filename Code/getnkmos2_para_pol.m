function [ncom,kcom]=getnkmos2_para_pol(wl)
% Returns the real and imaginary part of refractive index of MoS2 when
% the polarisation is parallel to its plane
% wl is in nm.
load wl01_nmos2_com_ponly_kmos2_com_refrac.mat
[r,~]=size(wl);
[~,loc]=ismember(wl,wl01);

ncom=nmos2_com(loc);
kcom=kmos2_com(loc);

if r==1
    ncom=ncom';
    kcom=kcom';
end
end