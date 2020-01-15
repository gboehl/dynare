% Stochastic growth model of Brock and Mirman (1972) with technology shock
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
% =========================================================================

var
    C        ${C}$ (long_name='consumption')
    K        ${K}$ (long_name='capital')
    A        ${Z}$ (long_name='total factor productivity')
;
varobs C;

varexo
    eps_A    ${\varepsilon_A}$ (long_name='TFP shock')
;
    
parameters
    alph    ${\alpha}$   (long_name='capital share')
    betta   ${\beta}$    (long_name='discount factor')
    rhoA    ${\rho_A}$   (long_name='persistence TFP')
    sigA    ${\sigma_A}$ (long_name='standard deviation TFP shock')
;

alph  = 0.35;
betta = 0.99;
rhoA  = 0.9;
sigA  = 0.6;

model;
[name='Euler equation']
C^(-1)=alph*betta*C(+1)^(-1)*A(+1)*K^(alph-1);
[name='capital law of motion'] 
K=A*K(-1)^alph-C;
[name='exogenous TFP process']
log(A)=rhoA*log(A(-1))+sigA*eps_A;
end;

shocks;
    var eps_A = 1;
end;

steady_state_model;
    A = 1;                           % technology level
    K = (alph*betta*A)^(1/(1-alph)); % capital level
    C = A*K^alph-K;                  % consumption level
end;

steady; % compute steady state given the starting values
resid;  % check the starting values for the steady state
check;  % check Blanchard & Khan rank condition

@#ifdef kronflag
identification(ar=3, useautocorr=1, nodisplay, nograph, parameter_set=calibration, analytic_derivation_mode=@{kronflag});
@#else
identification(ar=3, useautocorr=1, nodisplay, nograph, parameter_set=calibration, analytic_derivation_mode=0);
@#endif
