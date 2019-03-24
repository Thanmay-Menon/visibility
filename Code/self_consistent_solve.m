function [np,nowtheta]=self_consistent_solve(befn,para_n,perp_n,beftheta)
% This function solves for the effective refractive index in a
% self-consistent manner as described in the paper. befn is the refractive
% index in the previous layer, para_n is the (usually complex, to be given in the form, n-ik)
% refractive index when the polarisation is parallel to its plane and perp_n 
% is the refractive index when the polarisation is perpendicular to its
% plane. beftheta is the angle of incidence to this layer or the angle of refraction in
% the previous layer (value of both angles will be same due to geometry)
% Returns the effective refractive index and the angle of refraction in the
% current layer. All angles are in radians
[r,c]=size(beftheta);
x0=zeros(r,c);
func=@(x)(temp_func(x,para_n,perp_n,beftheta,befn));
options = optimoptions('fsolve','Display','none');
[nowtheta,~,e]=fsolve(func,x0,options);
if e==1
    np=1/sqrt(cos(nowtheta).^2./para_n.^2+sin(nowtheta).^2./perp_n.^2);
else
    display("Problem in fsolve in self_consistent_solve.m");
end
end

function out=temp_func(x,para_n,perp_n,beftheta,befn)
n_trial=1/sqrt(cos(x).^2./para_n.^2+sin(x).^2./perp_n.^2);

out=x-asin(befn.*sin(beftheta)./n_trial);
end