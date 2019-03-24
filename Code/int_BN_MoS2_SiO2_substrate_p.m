function iout=int_BN_MoS2_SiO2_substrate_p(wl,d_grap,d_bn,d_mos2,d_sio2,angle)

d_grap=d_grap.*10^-9;
d_bn=d_bn.*10^-9;
d_mos2=d_mos2.*10^-9;
d_sio2=d_sio2.*10^-9;

[Re_si,Im_si]=getnsi(wl);
Im_gra=getkgra(wl);
[nmos2_para,kmos2_para]=getnkmos2_para_pol(wl);
nmos2_perp=getnmos2_perp_pol(wl);

n0=1;
n1=1;
n2=2.2;
n4=1.47;
n5=Re_si-1i.*Im_si;

d1=d_grap;
d2=d_bn;
d3=d_mos2;
d4=d_sio2;

angle0=angle;
angle1=angle0;
angle2=get_angle(n1,n2,angle1);
[n3,angle3]=self_consistent_solve(n2,nmos2_para-1i.*kmos2_para,nmos2_perp,angle2);
angle4=get_angle(n3,n4,angle3);
angle5=get_angle(n4,n5,angle4);

wl=wl.*10^-9;

p1=4.*pi.*n1.*cos(angle1).*d1./wl;
p2=4.*pi.*n2.*cos(angle2).*d2./wl;
p3=4.*pi.*n3.*cos(angle3).*d3./wl;
p4=4.*pi.*n4.*cos(angle4).*d4./wl;

r1=(n1.*cos(angle0)-n0.*cos(angle1))./(n0.*cos(angle1)+n1.*cos(angle0));
r2=(n2.*cos(angle1)-n1.*cos(angle2))./(n1.*cos(angle2)+n2.*cos(angle1)); 
r3=(n3.*cos(angle2)-n2.*cos(angle3))./(n2.*cos(angle3)+n3.*cos(angle2));
r4=(n4.*cos(angle3)-n3.*cos(angle4))./(n3.*cos(angle4)+n4.*cos(angle3));
r5=(n5.*cos(angle4)-n4.*cos(angle5))./(n4.*cos(angle5)+n5.*cos(angle4));
    
iout=abs((r5+r4.*exp(1i.*p4)+r3.*exp(1i.*(p3+p4))+r2.*exp(1i.*(p2+p3+p4))+r1.*r2.*r5.*exp(1i.*p1)+r2.*r3.*r5.*exp(1i.*p2)+r3.*r4.*r5.*exp(1i.*p3)+r1.*r2.*r4.*exp(1i.*(p1+p4))+r1.*r3.*r5.*exp(1i.*(p1+p2))+r2.*r3.*r4.*exp(1i.*(p2+p4))+r2.*r4.*r5.*exp(1i.*(p2+p3))+r1.*exp(1i.*(p1+p2+p3+p4))+r1.*r2.*r3.*exp(1i.*(p1+p3+p4))+r1.*r3.*r4.*exp(1i.*(p1+p2+p4))+r1.*r4.*r5.*exp(1i.*(p1+p2+p3))+r1.*r2.*r3.*r4.*r5.*exp(1i.*(p1+p3)))./ ...
(r1.*r5+r2.*r5.*exp(1i.*p1)+r1.*r4.* exp(1i.*p4)+exp(1i.*(p1+p2+p3+p4))+r1.*r3.*exp(1i.*(p3+p4))+r2.*r4.*exp(1i.*(p1+p4))+r3.*r5.*exp(1i.*(p1+p2))+r1.*r2.*r3.*r5.*exp(1i.*(p2))+r1.*r3.*r4.*r5.*exp(1i.*p3)+r1.*r2.*exp(1i.*(p2+p3+p4))+r2.*r3.*exp(1i.*(p1+p3+p4))+r3.*r4.*exp(1i.*(p1+p2+p4))+r4.*r5.*exp(1i.*(p1+p2+p3))+r1.*r2.*r3.*r4.*exp(1i.*(p2+p4))+r1.*r2.*r4.*r5.*exp(1i.*(p2+p3))+r2.*r3.*r4.*r5.*exp(1i.*(p1+p3)))).^2;
end

function angle_curr=get_angle(n_pre,n_curr,angle_pre)
angle_curr=asin(n_pre.*sin(angle_pre)./n_curr);
end