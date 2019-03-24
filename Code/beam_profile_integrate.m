function out=beam_profile_integrate(intensity,NA,varargin)
%Assumes gaussian beam profile and integrates the intensity over angles
%from 0 to asin(NA) and returns the result of this integration.
% Takes precision (in radians) of the integration variable as an optional argument 
% Default precision is 0.001 radians
%intensity is in the form of a row vector and for the code to run without
%errors, length(intensity) has to be equal to length(ang_vec) (See Below)
prec=0.001;
if length(varargin)>1
    display('ERROR. At most 1 optional argument allowed');
elseif length(varargin)==1
    prec=varargin{1};
end
ang_vec=0:prec:round(asin(NA),3);% Would throw up an error if round(asin(NA),3) is not a multiple of prec
out=trapz(ang_vec,tan(ang_vec).*gaussian_beam(ang_vec,NA).*intensity);%trapezoidal integration is performed here
end

function beam_out=gaussian_beam(angle,NA)
dev=asin(NA);
beam_out=exp(-2.*(tan(angle)).^2./(tan(dev).^2));
end