% =========================================================================
% Stochastic growth model of Brock and Mirman (1972) with technology shock
% Willi Mutschler, January 2018
% willi@mutschler.eu
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