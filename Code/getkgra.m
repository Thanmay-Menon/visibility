function kg=getkgra(wl)
%Returns the imaginary part of the refractive index of graphene when the
%polarisation is parallel to its plane
% wl is in nm
load kgra_wl01

[r,~]=size(wl);
[~,loc]=ismember(wl,wl01);
kg=kgra(loc);
if r==1
    kg=kg';
end
end