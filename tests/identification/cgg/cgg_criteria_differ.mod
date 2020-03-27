% This is the Clarida, Gali and Gertler Basic New Keynesian model
% This mod file illustrates that due to numerical errors and numerical
% settings the identification criteria might differ
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
y        ${y}$             (long_name='output')
c        ${c}$             (long_name='consumption')
lam      ${\lambda}$       (long_name='marginal utility, i.e. lagrange multiplier budget')  
R        ${R}$             (long_name='nominal interest rate')  
pie      ${\pi}$           (long_name='inflation rate')    
n        ${n}$             (long_name='labor')
w        ${w}$             (long_name='real wage')
mc       ${mc}$            (long_name='marginal costs')
pstar    ${p^\ast}$        (long_name='price dispersion')
ptilde   ${\tilde{p}}$     (long_name='reoptimized price')
S1       ${S_{1}}$         (long_name='auxiliary variable for nonlinear pricing 1')
S2       ${S_{2}}$         (long_name='auxiliary variable for nonlinear pricing 2')
a        ${a}$             (long_name='total factor productivity')
RR       ${R^{**}}$        (long_name='natural interest rate')
yy       ${y^{**}}$        (long_name='natural output')
x        ${x}$             (long_name='output gap')
da       ${\delta a}$      (long_name='technology growth')
;

varexo
epsa     ${\varepsilon^A}$      (long_name='TFP shock')
epsR     ${\varepsilon^R}$      (long_name='monetary policy shock')
;

parameters
BETA     ${\beta}$         (long_name='discount factor')
PIESTAR  ${\pi^\ast}$      (long_name='annual target inflation rate')
PSIPIE   ${\psi_\pi}$      (long_name='Taylor rule parameter inflation')
PSIY     ${\psi_y}$        (long_name='Taylor rule parameter output')  
RHOR     ${\rho_R}$        (long_name='persistence Taylor rule')
SIGR     ${\sigma_R}$      (long_name='standard deviation monetary policy shock')
ETA      ${\eta}$          (long_name='elasticity of substitution')
THETA    ${\theta}$        (long_name='Calvo probability of resetting prices')
PSI      ${\psi}$          (long_name='labor disutility parameter')
RHOA     ${\rho_A}$        (long_name='persistence TFP')
SIGA     ${\sigma_A}$      (long_name='standard deviation TFP shock')
SUBSIDY  ${\nu}$           (long_name='subsidy')
;

%Calibration
BETA     = 0.99;
PIESTAR  = 1;
PSIPIE   = 1.5;
PSIY     = 0;
RHOR     = 0;
SIGR     = 0.2;
ETA      = 2;
THETA    = 3/4;
PSI      = 1;
RHOA     = 0.5;
SIGA     = 0.01;
SUBSIDY  = 1 - (ETA-1)/ETA;

model;

[name='marginal utility of consumption']
lam = 1/c;

[name='Euler equation']
lam = BETA*R/pie(+1)*lam(+1);

[name='labor supply']
lam*w = n^PSI;

[name='aggregate supply']
pstar*y = a*n;

[name='labor demand']
w*(1-SUBSIDY) = mc*a;

[name='market clearing']
y = c;

[name='recursive Price']
ptilde*S1 = ETA/(ETA-1)*S2;

[name='recursive Price 1']
S1 = y + THETA*BETA*lam(+1)/lam*pie(+1)^(ETA-1)*S1(+1);

[name='recursive Price 2']
S2 = y*mc + THETA*BETA*lam(+1)/lam*pie(+1)^(ETA)*S2(+1);

[name='aggregate price index']
1 = (1-THETA)*ptilde^(1-ETA) + THETA*pie^(ETA-1);

[name='overall price dispersion']
pstar = (1-THETA)*ptilde^(-ETA) + THETA*pie^(ETA)*pstar(-1);

[name='Taylor rule']
R = steady_state(R)^(1-RHOR)*R(-1)^RHOR*(pie/PIESTAR)^PSIPIE*(y/steady_state(y))^PSIY*exp(epsR);

[name='evolution of technology']
log(a) = RHOA*log(a(-1)) + epsa;

[name='technology growth']
da = log(a/a(-1));

[name='efficient interest rate']
RR = 1/BETA * a(+1)/a;

[name='efficient output']
yy = a;

[name='output gap']
x = (y-yy);

end;

steady_state_model;
a = 1;
da = 0;
pie = PIESTAR;
R = pie/BETA;
ptilde = ( (1-THETA*pie^(ETA-1))/ (1-THETA) )^(1/(1-ETA));
pstar = ( (1-THETA) / (1-THETA*pie^(ETA)) )*ptilde^(-ETA);
mc = (ETA-1)/ETA*ptilde* (1-THETA*BETA*pie^(ETA)) / (1-THETA*BETA*pie^(ETA-1));
w = 1/(1-SUBSIDY)*mc*a;
n = (pstar*mc/(1-SUBSIDY))^(1/(PSI+1));
y = pstar^(-1)*a*n;
c = y;
lam = 1/c;
S1 = y/(1-THETA*BETA*pie^(ETA-1));
S2 = y*mc/(1-THETA*BETA*pie^(ETA));
RR = 1/BETA;
yy = a;
x = y-yy;
end;


shocks;
var epsa = SIGA^2;
var epsR = SIGR^2;
end;


model_diagnostics;
steady;
resid;
check;

estimated_params;
PIESTAR  , 1;
PSIPIE   , 1.5;
PSIY     , 0;
RHOR     , 0;
ETA      , 1.5;
end;

varobs y c;
identification(tol_rank=1e-13, tol_sv=1e-2, tol_deriv=1e-6, normalize_jacobians=0, checks_via_subsets=0);
identification(tol_rank=1e-13, tol_sv=1e-2, tol_deriv=1e-6, normalize_jacobians=1, checks_via_subsets=0);
identification(tol_rank=1e-13, tol_sv=1e-2, tol_deriv=1e-6, normalize_jacobians=0, checks_via_subsets=1);
identification(tol_rank=1e-13, tol_sv=1e-2, tol_deriv=1e-6, normalize_jacobians=1, checks_via_subsets=1);

estim_params_ = [];
identification(checks_via_subsets=1, max_dim_subsets_groups=4);
