% Created by Marco Ratto (@rattoma, marco.ratto@ec.europa.eu)
% based on Kim, Jinill, 2003. "Functional equivalence between intertemporal and 
% multisectoral investment adjustment costs," Journal of Economic Dynamics 
% and Control, 27(4), pages 533-549.
% =========================================================================
% Copyright © 2010-2020 Dynare Team
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

var c k i a lam;
varexo ea;

parameters
alph
betae
delta
as
s
theta
phi
dumpy
;


alph=0.6;
betae=0.99;
delta=0.0125;
as=0.1;
phi=2;
s=betae*delta*alph/(1-betae+delta*betae);
theta=1;

model;
//(1-(betae*delta*alph/(1-betae+delta*betae)))^theta/c^(1+theta)/(1+theta)*(1+theta)*(i/(betae*delta*alph/(1-betae+delta*betae)))^theta*(i/k/delta)^phi=betae*(1-(betae*delta*alph/(1-betae+delta*betae)))^theta/c(+1)^(1+theta)/(1+theta)*(alph*(1+theta)*a(+1)^(1+theta)*k^(alph*(1+theta)-1)+(1-delta)*(i(+1)/k/delta)^phi*(1+theta)*(i(+1)/(betae*delta*alph/(1-betae+delta*betae)))^theta);
lam*(1+theta)*(i/(betae*delta*alph/(1-betae+delta*betae)))^theta*(i/k/delta)^phi=betae*lam(+1)*(alph*(1+theta)*a(+1)^(1+theta)*k^(alph*(1+theta)-1)+(1-delta)*(i(+1)/k/delta)^phi*(1+theta)*(i(+1)/(betae*delta*alph/(1-betae+delta*betae)))^theta);
k=(delta*(i/delta)^(1-phi)+(1-delta)*k(-1)^(1-phi))^(1/(1-phi));
((1-(betae*delta*alph/(1-betae+delta*betae)))*(c/(1-(betae*delta*alph/(1-betae+delta*betae))))^(1+theta) + 
(betae*delta*alph/(1-betae+delta*betae))*(i/(betae*delta*alph/(1-betae+delta*betae)))^(1+theta))^(1/(1+theta))=
(a*k(-1)^alph);
a = as+ea;
lam = (1-(betae*delta*alph/(1-betae+delta*betae)))^theta/c^(1+theta)/(1+theta);
//cobs = c+ec;
end;

steady_state_model;
s=betae*delta*alph/(1-betae+delta*betae);
a=as; %as^((1-alph)/(1+theta))*(delta^((phi+theta+1)/(theta+1))/s)^alph;
k=(delta/s/a)^(1/(alph-1));
i=delta*k;
c=(((a*k^alph)^(1+theta)-s*(i/s)^(1+theta))/(1-s))^(1/(1+theta))*(1-s);
lam = (1-s)^theta/c^(1+theta)/(1+theta);
end;

steady;
check;

shocks;
var ea = 1;
//var ec = 0;
end;

estimated_params;
alph ,uniform_pdf, , ,0.5,0.7;
//betae ,uniform_pdf,0.99,0.004,0.98,1;
//delta ,uniform_pdf,0.0125,0.001,0.01,0.015;
phi ,uniform_pdf, , ,0,10;
theta ,uniform_pdf, , ,0,10;
dumpy ,uniform_pdf, , ,0,10;
end;

varobs c i;

identification(advanced=1,max_dim_cova_group=3,tol_rank=1e-8);
//varobs c i lam; //to check if observing lam identifies phi and theta
//identification(ar=1,advanced=1,max_dim_cova_group=3,prior_mc=250);
//identification(prior_mc=100);

% Unit test for analytic_derivation_mode
load('kim2/identification/kim2_prior_mean_identif.mat','store_options_ident')
if store_options_ident.analytic_derivation~=1 && store_options_ident.analytic_derivation_mode~=-2
    error('the steady state file changed parameters and we should switch to numerical derivatives for the steady state, i.e. analytic_derivation_mode=-2')
end

% Unit test for correct identification results
load('kim2/identification/kim2_prior_mean_identif.mat','ide_moments_point', 'ide_spectrum_point', 'ide_minimal_point', 'ide_reducedform_point')
pause(1);
chk.ind0       = [1 1 1 0];
chk.indno      = [0 0 0 1; 0 1 1 0];
chk.jweak      = [0 1 1 0];
chk.jweak_pair = [0 0 0 0 1 0 0 0 0 0];
for strVars = {'ind0' 'indno' 'jweak' 'jweak_pair'}
    if ~isequal(ide_moments_point.(strVars{:}) , chk.(strVars{:}))
        disp('dMoments:')
        disp(ide_moments_point.dMOMENTS);
        disp(strVars{:})
        disp(ide_moments_point.(strVars{:}));
        error('identification based on moments is wrong for %s',strVars{:})
    end
    if ~isequal(ide_spectrum_point.(strVars{:}) , chk.(strVars{:}))
        disp('dSPECTRUM');
        disp(ide_spectrum_point.dSPECTRUM);
        disp(strVars{:})
        disp(ide_spectrum_point.(strVars{:}));
        error('identification based on spectrum is wrong for %s',strVars{:})
    end
    if ~isequal(ide_minimal_point.(strVars{:}) , chk.(strVars{:}))
        disp('dMINIMAL')
        disp(ide_minimal_point.dMINIMAL);
        disp(strVars{:})
        disp(ide_minimal_point.(strVars{:}));
        error('identification based on minimal system is wrong for %s',strVars{:})
    end
end

% Integration test if identification works without priors
estim_params_=[]; 
dumpy=0;
identification(advanced=1,max_dim_cova_group=3);
