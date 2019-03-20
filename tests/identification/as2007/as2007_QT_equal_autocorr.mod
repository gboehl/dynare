% this is the exact model Qu and Tkachenk (2012, Quantitative Economics)
% used in their replication file.
% This file illustrates that identification criteria are only local conditions,
% as when setting all autocorrelation coefficients to 0.5, 
% [psi1,psi2] are not identified instead of [psi1,psi2,rhor,sig2r].
% Created by Willi Mutschler (willi@mutschler.eu)
var z g R y pie c piep yp;
varexo e_z e_g e_R ;
parameters tau betta nu phi pistar psi1 psi2 rhor rhog rhoz sig2r sig2g sig2z;

tau = 2.0000;
betta = 0.9975;
nu = 0.1000;
phi= 53.6797;
pistar=1.0080^2; %note pistar is squared
psi1 = 1.5000;
psi2 = 0.1250;
rhor = 0.500;
rhog = 0.500;
rhoz = 0.500;
sig2r= (.002^2)*(10^5);
sig2g= (.006^2)*(10^5);
sig2z=(.003^2)*(10^5);


model;
#kap=(tau*(1-nu))/(nu*phi*pistar);
y = y(+1) + g - g(+1) - 1/tau*(R-pie(+1)-z(+1));
pie = betta*pie(+1)+kap*(y-g);
c = y - g;
R = rhor*R(-1)+(1-rhor)*psi1*pie+(1-rhor)*psi2*(y-g)+sqrt(sig2r)*e_R;
g = rhog*g(-1)+sqrt(sig2g)*e_g;
z = rhoz*z(-1)+sqrt(sig2z)*e_z;
piep = pie(+1);
yp = y(+1);
end;

estimated_params;
tau , 2.0000;
betta , 0.9975;
nu , 0.1000;
phi , 53.6797;
pistar ,1.0080^2;
psi1 , 1.5000;
psi2 , 0.1250;
rhor , 0.500;
rhog , 0.500;
rhoz , 0.500;
sig2r , (.002^2)*(10^5);
sig2g , (.006^2)*(10^5);
sig2z , (.003^2)*(10^5);
end;

steady_state_model;
pie=0; y=0; c=0; R=0; g=0; z=0; yp=0; piep=0;
end;
steady;
check;

varobs R y pie c;

shocks;
var e_z;
stderr 1;
var e_g;
stderr 1;
var e_R;
stderr 1;
end;

identification;
identification(useautocorr=0,ar=10);