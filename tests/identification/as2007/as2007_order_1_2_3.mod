%this is the mod file used in replication files of An and Schorfheide (2007)
% modified to include some obvious and artificial identification failures
% and to check whether all kronflags are working
% created by Willi Mutschler (@wmutschl, willi@mutschler.eu)
% =========================================================================
% Copyright (C) 2019-2020 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
% =========================================================================

var y R g z c dy p YGR INFL INT;
varobs y R g z c dy p YGR INFL INT;
varexo e_r e_g e_z;
parameters sigr sigg sigz tau phi psi1 psi2 rhor rhog rhoz rrst pist gamst nu cyst dumpy dumpyrhog;

rrst = 1.0000;
pist = 3.2000;
gamst= 0.5500;
tau  = 2.0000; 
nu   = 0.1000;
kap  = 0.3300;
phi  = tau*(1-nu)/nu/kap/exp(pist/400)^2;
cyst = 0.8500;
psi1 = 1.5000;
psi2 = 0.1250;
rhor = 0.7500;
rhog = 0.9500;
rhoz = 0.9000;
sigr = 0.2;
sigg = 0.6;
sigz = 0.3;
dumpy = 0;
dumpyrhog = 1;

model;
#pist2 = exp(pist/400);
#rrst2 = exp(rrst/400);
#bet   = 1/rrst2;
#gst   = 1/cyst;
#cst   = (1-nu)^(1/tau);
#yst   = cst*gst;
1 = exp(-tau*c(+1)+tau*c+R-z(+1)-p(+1));
(1-nu)/nu/phi/(pist2^2)*(exp(tau*c)-1) = (exp(p)-1)*((1-1/2/nu)*exp(p)+1/2/nu) - bet*(exp(p(+1))-1)*exp(-tau*c(+1)+tau*c+dy(+1)+p(+1));
exp(c-y) = exp(-g) - phi*pist2^2*gst/2*(exp(p)-1)^2;
R = rhor*R(-1) + (1-rhor)*psi1*p + (1-rhor)*psi2*(y-g) + sigr*e_r;
g = dumpyrhog*rhog*g(-1) + sigg*e_g;
z = rhoz*z(-1) + sigz*e_z;
YGR = gamst+100*(dy+z);
INFL = pist+400*p;
INT = pist+rrst+4*gamst+400*R;
dy = y - y(-1);
end;

shocks;
var e_r = 0.6^2;
var e_g = 0.5^2;
var e_z = 0.4^2;
corr e_r, e_g = 0.3;
corr e_r, e_z = 0.2;
corr e_z, e_g = 0.1;
end;

steady_state_model;
z=0; g=0; c=0; y=0; p=0; R=0; dy=0;
YGR=gamst; INFL=pist; INT=pist+rrst+4*gamst;
end;

estimated_params;
tau,   2,     1e-5, 10,      gamma_pdf,     2,    0.5; 

%these parameters do not enter the linearized solution
cyst, 0.85,  1e-5, 0.99999, beta_pdf,      0.85, 0.1;
sigg, 0.6,   1e-8, 5,     inv_gamma_pdf, 0.4,  4;
rhoz,  0.9,   1e-5, 0.99999, beta_pdf,      0.66, 0.15;
corr e_r,e_g, 0.3,   1e-8, 5, inv_gamma_pdf, 0.4,  4;
corr e_z,e_g, 0.3,   1e-8, 5, inv_gamma_pdf, 0.4,  4;
corr e_z,e_r, 0.3,   1e-8, 5, inv_gamma_pdf, 0.4,  4;

%these parameters could only be identified from the steady state of YGR INFL and INT, however, we observer y pi R instead
rrst,  1,     1e-5, 10,      gamma_pdf,     0.8,  0.5;
gamst, 0.55,  -5,   5,       normal_pdf,    0.4,  0.2;
dumpy, 0, -10, 10, normal_pdf, 0, 1;

%these parameters jointly determine the slope kappa of the linearized new keynesian phillips curve
pist,  3.2,   1e-5, 20,      gamma_pdf,     4,    2;
nu,    0.1,   1e-5, 0.99999, beta_pdf,      0.1,  .05;
phi,   50,    1e-5, 100,     gamma_pdf,     50,   20;

%these parameters are pairwise collinear as one should not use both formulations for the standard error of a shock
sigz, 0.3,   1e-8, 5,     inv_gamma_pdf, 0.4,  4;
stderr e_z, 0.3,   1e-8, 5, inv_gamma_pdf, 0.4,  4;

%these parameters are pairwise collinear as they are multiplicative
rhog,  0.95,  1e-5, 0.99999, beta_pdf,      0.8,  0.1;
dumpyrhog, 1, -10, 10, normal_pdf, 1, 1;

%these parameters are jointly not identified due to the specification of the Taylor rule
psi1,  1.5,   1e-5, 10,      gamma_pdf,     1.5,  0.25;
psi2,  0.125, 1e-5, 10,      gamma_pdf,     0.5,  0.25;
rhor,  0.75,  1e-5, 0.99999, beta_pdf,      0.5,  0.2;
stderr e_r, 0.2,   1e-8, 5, inv_gamma_pdf, 0.3,  4;

end;

steady;
check;
stoch_simul(order=3,irf=0); %needed for identification(order=3)

@#for ORDER in [1, 2, 3]
@#for KRONFLAG in [-1, -2, 0]
fprintf('*** ORDER = @{ORDER} WITH ANALYTIC_DERIVATION_MODE=@{KRONFLAG} ***\n')
identification(order=@{ORDER}, parameter_set=calibration, grid_nbr=10,analytic_derivation_mode=@{KRONFLAG});
@#endfor
@#endfor