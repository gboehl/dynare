// Uses autonomous system from sim_base.mod, but adds separate system where exogenous variables have several leads and lags
// Lags and leads on exogenous variables are substituted out by auxiliary variables

var c k z_backward z_forward x_lag_1 x_lag_2 x_lag_3 x_lag_4 x_lead_1 x_lead_2 x_lead_3 x_lead_4;
varexo x;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;

model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1); // Resource constraint
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam); // Euler equation
z_backward=0.1*1+0.9*z_backward(-1)+(x_lag_4-1);
z_forward=0.2*1+0.8*z_forward(+1)+(x_lead_4-1);

x_lag_1=x(-1);
x_lag_2=x_lag_1(-1);
x_lag_3=x_lag_2(-1);
x_lag_4=x_lag_3(-1);
x_lead_1=x(+1);
x_lead_2=x_lead_1(+1);
x_lead_3=x_lead_2(+1);
x_lead_4=x_lead_3(+1);
end;

initval;
c = 1.2;
k = 12;
x = 1; %set x(0)
x_lag_1 = 1.3; %sets x(-1)
x_lag_2 = 1.3; %sets x(-2)
end;

shocks;
var x; %sets x(+2)
periods 2;
values 0.9;
end;

simul(periods=200,maxit=100);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

base_results=load('sim_base_results.mat');
if max(abs(base_results.oo_.endo_simul(strmatch('c',base_results.M_.endo_names,'exact'),1+base_results.M_.maximum_endo_lag:end-base_results.M_.maximum_endo_lead) -...
    oo_.endo_simul(strmatch('c',M_.endo_names,'exact'),1+M_.maximum_endo_lag:end-M_.maximum_endo_lead)))>1e-8 || ...
    max(abs(base_results.oo_.endo_simul(strmatch('k',base_results.M_.endo_names,'exact'),1+base_results.M_.maximum_endo_lag:end-base_results.M_.maximum_endo_lead) -...
        oo_.endo_simul(strmatch('k',M_.endo_names,'exact'),1+M_.maximum_endo_lag:end-M_.maximum_endo_lead)))>1e-8
    error('Autonomous system part is wrong')
end

