function iout=new_stuff(wl,d_grap,d_bn,d_sio2,angle,pol)
d_grap=d_grap.*10^-9;
d_bn=d_bn.*10^-9;
d_sio2=d_sio2.*10^-9;
%n0 is air,n1 is graphene,thixb nckness 0.34nm,n2 is BN,n3 is SiO2,n4 is Si
[Re_si,Im_si]=getnsi(wl);
Im_gra=getkgra(wl);

n0=1;
n1=2.7-1i.*Im_gra;
n2=2.2;
n3=1.47;
n4=Re_si-1i.*Im_si;

d1=d_grap;
d2=d_bn;
d3=d_sio2;

angle0=angle;
angle1=get_angle(1,n1,angle0);
angle2=get_angle(n1,n2,angle1);
angle3=get_angle(n2,n3,angle2);
angle4=get_angle(n3,n4,angle3);

wl=wl.*10^-9;



p1=4.*pi.*n1.*cos(angle1).*d1./wl;
p2=4.*pi.*n2.*cos(angle2).*d2./wl;
p3=4.*pi.*n3.*cos(angle3).*d3./wl;

if(pol=='s')
    r1=(n0.*cos(angle0)-n1.*cos(angle1))./(n0.*cos(angle0)+n1.*cos(angle1));
    r2=(n1.*cos(angle1)-n2.*cos(angle2))./(n1.*cos(angle1)+n2.*cos(angle2));
    r3=(n2.*cos(angle2)-n3.*cos(angle3))./(n2.*cos(angle2)+n3.*cos(angle3));
    r4=(n3.*cos(angle3)-n4.*cos(angle4))./(n3.*cos(angle3)+n4.*cos(angle4));

elseif(pol=='p')
    r1=(n0.*cos(angle1)-n1.*cos(angle0))./(n0.*cos(angle1)+n1.*cos(angle0));
    r2=(n1.*cos(angle2)-n2.*cos(angle1))./(n1.*cos(angle2)+n2.*cos(angle1)); 
    r3=(n2.*cos(angle3)-n3.*cos(angle2))./(n2.*cos(angle3)+n3.*cos(angle2));
    r4=(n3.*cos(angle4)-n4.*cos(angle3))./(n3.*cos(angle4)+n4.*cos(angle3)); 
end

iout=(abs((r1+r2.*exp(-1i.*p1)+r3.*exp(-1i.*(p1+p2))+r4.*exp(-1i.*(p1+p2+p3))+r1.*r2.*r3.*exp(-1i.*p2)+r1.*r3.*r4.*exp(-1i.*p3)+r1.*r2.*r4.*exp(-1i.*(p3+p2))+r2.*r3.*r4.*exp(-1i.*(p1+p3)))./ ...
    (1+r1.*r2.*exp(-1i.*p1)+r1.*r3.*exp(-1i.*(p1+p2))+r1.*r4.*exp(-1i.*(p1+p2+p3))+r2.*r3.*exp(-1i.*p2)+r3.*r4.*exp(-1i.*p3)+r2.*r4.*exp(-1i.*(p2+p3))+r1.*r2.*r3.*r4.*exp(-1i.*(p1+p3))))).^2;
end

function angle_curr=get_angle(n_pre,n_curr,angle_pre)
angle_curr=asin(n_pre.*sin(angle_pre)./n_curr);
end