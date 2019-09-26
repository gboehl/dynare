% this is the exact model Qu and Tkachenk (2012, Quantitative Economics)
% used in their replication file.
% This file is used to check whether the G matrix is computed correctly.
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
rhor = 0.7500;
rhog = 0.9500;
rhoz = 0.9000;
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
tau,    uniform_pdf, , , 1.5, 2.5;
betta,  uniform_pdf, , , 0.95, 0.99999;
nu,     uniform_pdf, , , 0.05, 0.5;
phi,    uniform_pdf, , , 10, 80;
pistar, uniform_pdf, , , 1, 1.1;
psi1,   uniform_pdf, , , 1.05, 2.5;
psi2,   uniform_pdf, , , 0, 1.5;
rhor,   uniform_pdf, , , 0, 1;
rhog,   uniform_pdf, , , 0, 1;
rhoz,   uniform_pdf, , , 0, 1;
sig2r,  uniform_pdf, , , 0.01, 1;
sig2g,  uniform_pdf, , , 0.01, 1;
sig2z,  uniform_pdf, , , 0.01, 1;
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

identification(parameter_set=calibration,
               grid_nbr=10000,
               no_identification_strength,
               no_identification_reducedform,
               no_identification_moments,
               checks_via_subsets=1,
               max_dim_subsets_groups=4);

%     load('G_QT'); %note that this is computed using replication files of Qu and Tkachenko (2012)
%     temp = load([M_.dname filesep 'identification' filesep M_.fname '_identif']);
%     G_dynare = temp.ide_spectrum_point.G;
% 
%     % Compare signs
%     if ~isequal(sign(G_dynare),sign(G_QT))
%         error('signs of normalized G are note equal');
%     end
% 
%     % Compare normalized versions
%     tilda_G_dynare = temp.ide_spectrum_point.tilda_G;
%     ind_G_QT = (find(max(abs(G_QT'),[],1) > temp.store_options_ident.tol_deriv));
%     tilda_G_QT = zeros(size(G_QT));
%     delta_G_QT = sqrt(diag(G_QT(ind_G_QT,ind_G_QT)));
%     tilda_G_QT(ind_G_QT,ind_G_QT) = G_QT(ind_G_QT,ind_G_QT)./((delta_G_QT)*(delta_G_QT'));
%     if ~isequal(rank(tilda_G_QT,temp.store_options_ident.tol_rank),rank(tilda_G_dynare,temp.store_options_ident.tol_rank))
%         error('ranks are not the same for normalized version')
%     end
% 
%     max(max(abs(abs(tilda_G_dynare)-abs(tilda_G_QT))))
%     norm(tilda_G_dynare - tilda_G_QT)
