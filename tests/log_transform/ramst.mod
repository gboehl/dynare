// Test for var(log) in a deterministic context

var(log) c k;
varexo x;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;


model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
end;

initval;
x = 1;
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k;
end;

steady;

check;

shocks;
var x;
periods 1;
values 1.2;
end;

perfect_foresight_setup(periods=200);
perfect_foresight_solver;

idx_c = strmatch('c', M_.endo_names);
idx_log_c = strmatch('LOG_c', M_.endo_names);
if isempty(idx_log_c)
  error('Log-transformed variable not created')
end

if max(abs(oo_.endo_simul(idx_log_c, :) - log(oo_.endo_simul(idx_c, :)))) > 1e-7
  error('Transformation not correctly performed')
end

S = load('../simul/ramst/Output/ramst_results.mat');
if max(abs(oo_.endo_simul(idx_c, :) - S.oo_.endo_simul(idx_c, :))) > 1e-7
  error('Result differs from non-transformed model')
end
