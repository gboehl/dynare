% Original model by: J. Linde, M. Trabandt (2019: Should We Use Linearized Models to Calculate Fiscal Multipliers?
% Journal of Applied Econometrics, 2018, 33: 937-965. http://dx.doi.org/10.1002/jae.2641
% This version has some additional dynamics for capital and investment
% Created by Willi Mutschler (@wmutschl, willi@mutschler.eu)
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
% =========================================================================
% Declare endogenous variables
% =========================================================================

var 
    % Staggered-price economy
    c        ${c}$              (long_name='consumption')
    lam      ${\lambda}$        (long_name='lagrange multiplier budget')
    n        ${n}$              (long_name='labor')
    w        ${w}$              (long_name='real wage')
    r        ${r}$              (long_name='interest rate')
    pie      ${\pi}$            (long_name='inflation')
    y        ${y}$              (long_name='output')
    pstar    ${p^\ast}$         (long_name='price dispersion')
    vartheta ${\vartheta}$      (long_name='aggregate price index')
    mc       ${mc}$             (long_name='marginal costs')
    s        ${s}$              (long_name='auxiliary variable for nonlinear pricing 1')
    f        ${f}$              (long_name='auxiliary variable for nonlinear pricing 2')
    a        ${a}$              (long_name='auxiliary variable for nonlinear pricing 3')
    ptilde   ${\tilde{p}}$      (long_name='reoptimized price')
    delta1   ${\Delta_1}$       (long_name='price dispersion 1')
    delta2   ${\Delta_2}$       (long_name='price dispersion 2') 
    delta3   ${\Delta_3}$       (long_name='price dispersion 3')
    b        ${b}$              (long_name='bonds') 
    tau      ${\tau}$           (long_name='lump-sum tax') 
    g
    nu %consumption preference shock

    % Flex-price economy
    cpot     ${c^{pot}}$        (long_name='flex-price consumption')
    npot     ${n^{pot}}$        (long_name='flex-price labor')
    wpot     ${w^{pot}}$        (long_name='flex-price real wage')
    rrpot     ${r^{pot}}$        (long_name='flex-price interest rate')
    ypot     ${y^{pot}}$        (long_name='flex-price output')
    bpot     ${b^{pot}}$        (long_name='flex-price bonds') 
    taupot   ${\tau^{pot}}$     (long_name='flex-price lump-sum tax') 

    % Added variables for capital and investment
    k        ${k}$              (long_name='capital')
    rk       ${r^{K}}$          (long_name='rental rate on capital')
    iv       ${i}$              (long_name='investment')
    q        ${q}$              (long_name='Tobins Q')    
    kpot     ${k^{pot}}$        (long_name='flex-price capital')
    rkpot    ${r^{K,pot}}$      (long_name='flex-price rental rate on capital')
    ivpot    ${i^{pot}}$        (long_name='flex-price investment')    
;

% =========================================================================
% Declare observable variables
% =========================================================================
varobs c y;

% =========================================================================
% Declare exogenous variables
% =========================================================================
varexo
epsg
epsnu
    ;

% =========================================================================
% Declare parameters
% =========================================================================
parameters
    ALPHA  ${\alpha}$   (long_name='capital share')
    BETA   ${\beta}$    (long_name='discount rate')
    PIEBAR ${\bar{\pi}}$    (long_name='target inflation rate')
    IOTA  ${\iota}$  (long_name='indexation')
    CHI   ${\chi}$    (long_name='inverse frisch elasticity')
    THETAP ${\theta_p}$    (long_name='net markup')
    XIP    ${\xi_p}$    (long_name='Calvo probability')
    PSI    ${\psi}$    (long_name='Kimbal curvature')
    GAMMAPIE ${\gamma_\pi}$    (long_name='Taylor rule parameter inflation')
    GAMMAX ${\gamma_x}$    (long_name='Taylor rule parameter output')
    DELTA  ${\delta}$   (long_name='depreciation rate')
    QBAR   ${\bar{Q}}$  (long_name='steady state Q')
    NUBAR  ${\bar{\nu}}$  (long_name='steady state consumption preference shock')
    BBAR  %government debt to annualized output ratio
    GBAR_o_YBAR  % government consumption ratio
    TAUBAR
    TAUNBAR
    PHI
    RHONU
    RHOG
    SIGNU
    SIGG
;


% =========================================================================
% Calibration
% =========================================================================
BETA = 0.995;
PIEBAR = 1 + 0.005;
ALPHA = 0.3;
NUBAR = 0;%0.01;
BBAR = 0;%0.6;
GBAR_o_YBAR = 0;%0.2;
TAUBAR = 0;
PHI = 0.0101;
IOTA = 0.5;
CHI = 1/0.4;
XIP = 2/3;
PSI = -12.2;
THETAP = 0.1;
GAMMAPIE = 1.5;
GAMMAX = 0.125;
QBAR=1;
DELTA = 0;%0.025;
RHONU = 0.95;
RHOG = 0.95;
SIGG = 0.6/100;
SIGNU = 0.6/100;
TAUNBAR = 0;

% ========================================================================
% Model equations
% =========================================================================
model;
% -------------------------------------------------------------------------
% Auxiliary expressions
% -------------------------------------------------------------------------
#OMEGAP = (1+THETAP)*(1+PSI)/(1+PSI+THETAP*PSI);   %gross markup
#gammap = pie^IOTA*PIEBAR^(1-IOTA);               %indexation rule
%#gammap = PIEBAR;
% -------------------------------------------------------------------------
% Added equations for capital and investment
% -------------------------------------------------------------------------
[name='foc household wrt i']
lam = lam*q;
[name='foc household wrt k']
lam*q = BETA*lam(+1)*(rk(+1)+ (1-DELTA)*q(+1) );
[name='capital accumulation']
k = (1-DELTA)*k(-1) + iv;
[name='capital demand']
rk=ALPHA*mc*y/k(-1);

[name='flex price foc household wrt k']
(cpot-steady_state(cpot)*nu)^(-1) = BETA*(cpot(+1)-steady_state(cpot)*nu(+1))^(-1)*(rkpot(+1) + (1-DELTA) );
[name='flex price capital accumulation']
kpot = (1-DELTA)*kpot(-1) + ivpot;
[name='flex price capital demand']
rkpot=ALPHA/(1+THETAP)*(npot/kpot(-1))^(1-ALPHA);

% -------------------------------------------------------------------------
% Actual model equations
% -------------------------------------------------------------------------
[name='foc household wrt c (marginal utility)']
(c-steady_state(c)*nu)^(-1) = lam;

[name='foc household wrt l (leisure vs. labor)']
n^CHI = (1-TAUNBAR)*lam*w;

[name='foc household wrt b (Euler equation)']
lam = BETA*r/pie(+1)*lam(+1);

[name='aggregate resource constraint']
y = c + g + iv;

[name='production function']
y=pstar^(-1)*k(-1)^ALPHA*n^(1-ALPHA);

[name='Nonlinear pricing 1']
s = OMEGAP*lam*y*vartheta^(OMEGAP/(OMEGAP-1))*mc + BETA*XIP*(gammap/pie(+1))^(OMEGAP/(1-OMEGAP))*s(+1);

[name='Nonlinear pricing 2']
f = lam*y*vartheta^(OMEGAP/(OMEGAP-1)) + BETA*XIP*(gammap/pie(+1))^(1/(1-OMEGAP))*f(+1);
  
[name='Nonlinear pricing 3']
a = PSI*(OMEGAP-1)*lam*y + BETA*XIP*(gammap/pie(+1))*a(+1);
  
[name='Nonlinear pricing 4']
s = f*ptilde -a*ptilde^(1+OMEGAP/(OMEGAP-1));

[name='zero profit condition']
vartheta = 1 + PSI - PSI*delta2;

[name='aggregate price index']
vartheta = delta3;

[name='overall price dispersion']
pstar = 1/(1+PSI)*vartheta^(OMEGAP/(OMEGAP-1))*delta1^(OMEGAP/(1-OMEGAP)) + PSI/(1+PSI);

[name='price dispersion 1']
delta1^(OMEGAP/(1-OMEGAP)) = (1-XIP)*ptilde^(OMEGAP/(1-OMEGAP)) + XIP*(gammap/pie*delta1(-1))^(OMEGAP/(1-OMEGAP));

[name='price dispersion 2']
delta2 = (1-XIP)*ptilde+XIP*gammap/pie*delta2(-1);

[name='price dispersion 3']
delta3^(1/(1-OMEGAP)) = (1-XIP)*ptilde^(1/(1-OMEGAP)) + XIP*(gammap/pie*delta3(-1))^(1/(1-OMEGAP));

[name='marginal costs']
(1-ALPHA)*mc = w*k(-1)^(-ALPHA)*n^ALPHA;

[name='Taylor Rule']
r = steady_state(r)*(pie/PIEBAR)^GAMMAPIE*( (y/steady_state(y))/(ypot/steady_state(ypot)) )^GAMMAX;

[name='government budget']
b = r(-1)/pie*b(-1) + g/steady_state(y) -TAUNBAR*w*n/steady_state(y) - tau;

[name='fiscal rule']
tau = TAUBAR + PHI*(b(-1)-BBAR);

[name='flex price Euler equation']
(cpot-steady_state(cpot)*nu)^(-1) = BETA*rrpot*(cpot(+1)-steady_state(cpot)*nu(+1))^(-1);

[name='flex price foc household wrt l (leisure vs. labor)']
npot^CHI = (1-TAUNBAR)*(cpot-steady_state(cpot)*nu)^(-1)*wpot;

[name='flex price wage']
(1-ALPHA)/(1+THETAP)*kpot(-1)^ALPHA = wpot*npot^ALPHA;

[name='flex price aggregate resource constraint']
ypot = cpot + g + ivpot;

[name='flex price production function']
ypot=kpot(-1)^ALPHA*npot^(1-ALPHA);

[name='flex price government budget']
bpot = rrpot(-1)*bpot(-1) + g/steady_state(y) -TAUNBAR*wpot*npot/steady_state(y) - taupot;

[name='fiscal rule ']
taupot = TAUBAR + PHI*(bpot(-1)-BBAR);

g/steady_state(y) - GBAR_o_YBAR = RHOG*(g(-1)/steady_state(y) - GBAR_o_YBAR) + epsg;
nu - NUBAR = RHONU*(nu(-1) - NUBAR) + epsnu;
end;


% =========================================================================
% Steady state using a steady_state_model block
% =========================================================================

steady_state_model;
OMEGP = (1+THETAP)*(1+PSI)/(1+PSI+THETAP*PSI);
q = QBAR;
pie = PIEBAR;
GAMMAP = PIEBAR^IOTA*PIEBAR^(1-IOTA);
b = BBAR;
tau = TAUBAR;
taun = TAUNBAR;
nu = NUBAR;
r = pie/BETA;
RK = (1/BETA+DELTA-1)*q; % foc k
aux1 = ( (1-XIP)/(1-XIP*(GAMMAP/pie)^(1/(1-OMEGP)) ) )^(1-OMEGP);
aux2 = (1-XIP)/(1-XIP*(GAMMAP/pie));
ptilde = (1+PSI)/ ( aux1 + PSI*aux2 );
delta2 = aux2*ptilde;
delta3 = aux1*ptilde;
vartheta = delta3;
delta1 = ( (1-XIP)/(1-XIP*(GAMMAP/pie)^(OMEGP/(1-OMEGP) )) )^((1-OMEGP)/OMEGP)*ptilde;
pstar = 1/(1+PSI)*vartheta^(OMEGP/(OMEGP-1))*delta1^(OMEGP/(1-OMEGP)) + PSI/(1+PSI);
f_o_a = ( vartheta^(OMEGP/(OMEGP-1)) / (1-BETA*XIP*(GAMMAP/pie)^(1/(1-OMEGP))) ) / ( PSI*(OMEGP-1) /(1-BETA*XIP*GAMMAP/pie) );
s_o_a = ( OMEGP*vartheta^(OMEGP/(OMEGP-1)) / (1-BETA*XIP*(GAMMAP/pie)^(OMEGP/(1-OMEGP))) ) / ( PSI*(OMEGP-1) /(1-BETA*XIP*GAMMAP/pie) );
mc = (f_o_a*ptilde-ptilde^(1+OMEGP/(OMEGP-1)))*s_o_a^(-1);

k_o_n = (RK/(mc*ALPHA))^(1/(ALPHA-1));
w = (1-ALPHA)*mc*k_o_n^ALPHA; % labor demand
y_o_n = pstar^(-1)*k_o_n^ALPHA; % production function
iv_o_n = DELTA*k_o_n;
iv_o_y = iv_o_n / y_o_n;
g_o_y = GBAR_o_YBAR;
c_o_n = y_o_n - iv_o_n - g_o_y; % market clearing

% The labor level, closed-form solution for labor
n = (c_o_n^(-1)*w)^(1/(CHI+1)) ;
% Value of remaining variables
c = c_o_n*n;
rk = RK;


iv = iv_o_n *n;
k = k_o_n*n;
y = y_o_n*n;
lam = n^CHI/w;
g = g_o_y*y;
a = PSI*(OMEGP-1)*y*lam/(1-BETA*XIP*GAMMAP/pie);
f = f_o_a*a;
s = f*ptilde - a*ptilde^(1+OMEGP/(OMEGP-1));

%flex price economy
bpot = BBAR;
taupot = TAUBAR;
taun = TAUNBAR;

rrpot = 1/BETA;
RKpot = (1/BETA+DELTA-1)*QBAR; % foc k
rkpot = RKpot;

kpot_o_npot = (RKpot/(ALPHA/(1+THETAP)))^(1/(ALPHA-1));
wpot = (1-ALPHA)/(1+THETAP)*kpot_o_npot^ALPHA; % labor demand
ypot_o_npot = kpot_o_npot^ALPHA; % production function
ivpot_o_npot = DELTA*kpot_o_npot;
ivpot_o_ypot = ivpot_o_npot / ypot_o_npot;
cpot_o_npot = ypot_o_npot - ivpot_o_npot - g_o_y; % market clearing
% The labor level, closed-form solution for labor
npot = (cpot_o_npot^(-1)*wpot)^(1/(CHI+1)) ;
% Value of remaining variables
cpot = cpot_o_npot*npot;
ivpot = ivpot_o_npot *npot;
kpot = kpot_o_npot*npot;
ypot = ypot_o_npot*npot;

end;



% =========================================================================
% Declare settings for shocks
% =========================================================================
shocks;
var epsg = 1;
var epsnu = 1;
end;

estimated_params;
ALPHA,       0.3,     normal_pdf, 0.3 ,    0.1;
BETA,        0.995,   normal_pdf, 0.995,   0.1;
PIEBAR,      1.005,   normal_pdf, 1.005,   0.1;
IOTA,        0.5,     normal_pdf, 0.5,     0.1;
CHI,         1/0.4,   normal_pdf, 1/0.4,   0.1;
THETAP,      0.1,     normal_pdf, 0.1,     0.1;
XIP,         2/3,     normal_pdf, 2/3,     0.1;
PSI,         -12.2,   normal_pdf, -12.2,   0.1;
GAMMAPIE,    1.5,     normal_pdf, 1.5,     0.1;
GAMMAX,      0.125,   normal_pdf, 0.125,   0.1;
DELTA,       0,       normal_pdf, 0,       0.1;
QBAR,        1,       normal_pdf, 1,       0.1;
NUBAR,       0,       normal_pdf, 0,       0.1;
BBAR,        0,       normal_pdf, 0,       0.1;
GBAR_o_YBAR, 0,       normal_pdf, 0,       0.1; %commenting this solves the analytic_derivation_mode=-2|-1 problem
TAUBAR,      0,       normal_pdf, 0,       0.1;
TAUNBAR,     0,       normal_pdf, 0,       0.1;
PHI,         0.0101,  normal_pdf, 0.0101,  0.1;
RHONU,       0.95,    normal_pdf, 0.95,    0.1;
RHOG,        0.95,    normal_pdf, 0.95,    0.1;
SIGNU,       0.6/100, normal_pdf, 0.6/100, 0.1;
SIGG,        0.6/100, normal_pdf, 0.6/100, 0.1;
end;


% =========================================================================
% Computations
% =========================================================================
model_diagnostics;
steady;
resid;
check;

identification(order=1,no_identification_strength,analytic_derivation_mode= 0,ar=5); %works
%identification(no_identification_strength,analytic_derivation_mode= 1,ar=5); %works, this takes the longest due to Kronecker products
options_.dynatol.x = 1e-9;
identification(order=1,no_identification_strength,analytic_derivation_mode=-1,ar=5); %works, but tolerance is way too small
identification(order=1,no_identification_strength,analytic_derivation_mode=-2,ar=5); %works, but tolerance is way to small
options_.dynatol.x = 1e-5; %this is the default value
identification(order=1,no_identification_strength,analytic_derivation_mode=-1,ar=5); %works only if GBAR_o_YBAR is uncommented
identification(order=1,no_identification_strength,analytic_derivation_mode=-2,ar=5); %works only if GBAR_o_YBAR is uncommented
