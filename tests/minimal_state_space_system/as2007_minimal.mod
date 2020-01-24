% this is the Smets and Wouters (2007) model for which Komunjer and Ng (2011)
% derived the minimal state space system. In Dynare, however, we use more 
% powerful minreal function
% created by Willi Mutschler (@wmutschl, willi@mutschler.eu)
% =========================================================================
% Copyright (C) 2020 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
% =========================================================================

var y R g z c dy p YGR INFL INT;
varobs y R p c YGR INFL INT;
varexo e_r e_g e_z;
parameters tau phi psi1 psi2 rhor rhog rhoz rrst pist gamst nu cyst;

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
R = rhor*R(-1) + (1-rhor)*psi1*p + (1-rhor)*psi2*(y-g) + e_r;
g = rhog*g(-1) + e_g;
z = rhoz*z(-1) + e_z;
YGR = gamst+100*(dy+z);
INFL = pist+400*p;
INT = pist+rrst+4*gamst+400*R;
dy = y - y(-1);
end;

shocks;
var e_r; stderr 0.2/100;
var e_g; stderr 0.6/100;
var e_z; stderr 0.3/100;
end;

steady_state_model;
z=0; g=0; c=0; y=0; p=0; R=0; dy=0;
YGR=gamst; INFL=pist; INT=pist+rrst+4*gamst;
end;
stoch_simul(order=1,irf=0,periods=0);
options_.qz_criterium = 1;

indx = [M_.nstatic+(1:M_.nspred)]';
indy = 1:M_.endo_nbr';

SS.A = oo_.dr.ghx(indx,:);
SS.B = oo_.dr.ghu(indx,:);
SS.C = oo_.dr.ghx(indy,:);
SS.D = oo_.dr.ghu(indy,:);

[CheckCO,minnx,minSS] = get_minimal_state_representation(SS,0);

Sigmax_full = lyapunov_symm(SS.A, SS.B*M_.Sigma_e*SS.B', options_.lyapunov_fixed_point_tol, options_.qz_criterium, options_.lyapunov_complex_threshold, 1, options_.debug);
Sigmay_full = SS.C*Sigmax_full*SS.C' + SS.D*M_.Sigma_e*SS.D';

Sigmax_min = lyapunov_symm(minSS.A, minSS.B*M_.Sigma_e*minSS.B', options_.lyapunov_fixed_point_tol, options_.qz_criterium, options_.lyapunov_complex_threshold, 1, options_.debug);
Sigmay_min = minSS.C*Sigmax_min*minSS.C' + minSS.D*M_.Sigma_e*minSS.D';

([Sigmay_full(:) - Sigmay_min(:)]')
sqrt(([diag(Sigmay_full), diag(Sigmay_min)]'))
dx = norm( Sigmay_full - Sigmay_min, Inf);
if dx > 1e-12
    error('something wrong with minimal state space computations')
else
    fprintf('numerical error for moments computed from minimal state system is %d\n',dx)
end

