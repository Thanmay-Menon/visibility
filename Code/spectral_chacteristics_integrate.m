function [outR,outG,outB]=spectral_chacteristics_integrate(intensity,varargin)
%intensity is in the form of a row vector and for the code to run without
%errors, length(intensity) has to be equal to length(wl_vec) (See Below)
%prec can be set as an optional argument. Default value is 1 nm. In my
%example, precision can't be higher than 0.1 nm as that is that is the
%precision of the spectral sensitivity table of my camera
prec=1;%prec is in nm
if prec<0.1
    display('Set prec to greater than 0.1'); 
elseif length(varargin)>1
    display('ERROR. At most 1 optional argument allowed');
elseif length(varargin)==1
    prec=varargin;
end

wl_vec=400:prec:750;% Would throw up an error if 750 is not a multiple of prec
[red_spec,green_spec,blue_spec]=get_spec(wl_vec);
outR=trapz(wl_vec,bbody3100K(wl_vec).*red_spec.*intensity);%trapezoidal integration is performed
outG=trapz(wl_vec,bbody3100K(wl_vec).*green_spec.*intensity);
outB=trapz(wl_vec,bbody3100K(wl_vec).*blue_spec.*intensity);
end


function [red_out,green_out,blue_out]=get_spec(wl)
load red_spec.mat
load green_spec.mat
load blue_spec.mat
load wl01.mat

[r,~]=size(wl);
[~,loc]=ismember(wl,wl01);

red=red_spec(:,2);
green=green_spec(:,2);
blue=blue_spec(:,2);

red_out=red(loc);
green_out=green(loc);
blue_out=blue(loc);

if r==1
   red_out=red_out';
green_out=green_out';
blue_out=blue_out';
end
end

function bsource=bbody3100K(wl)
hc=(6.626*10^-34)*(3*10^8);
kT=(1.3806488*10^-23)*3100;
wl=wl*10^-9;
bsource=(2*hc*(3*10^8))./((wl.^5).*(exp(hc./(wl.*kT))-1));
end

