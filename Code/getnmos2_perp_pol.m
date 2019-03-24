function n_perp=getnmos2_perp_pol(wl)
%Returns the refractive index of MoS2 when polarisation is perpendicular to the plane. 
% wl is in nm.
load wl01_nmos2_com_ponly_kmos2_com_refrac.mat

[r,~]=size(wl);
[~,loc]=ismember(wl,wl01);
n_perp=nmos2_ponly(loc);

if r==1
 n_perp=n_perp';
end
end