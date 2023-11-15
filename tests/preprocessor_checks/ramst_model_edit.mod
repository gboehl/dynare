/* Test for multiple model blocks, model_remove, model_options and var_remove commands,
   and model_replace block.
   It should give the same results as ramst.mod. */

var c k;
varexo x;

var dummy1 dummy2 dummy3;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;

model;
  [ name = 'ressource constraint' ]
  c + k = aa*x*k(-1)^alph;
end;

model;
  [ name = 'eq:dummy1', endogenous = 'dummy1' ]
  dummy1 = c + 1;
  [ foo = 'eq:dummy2' ] // Since dummy2 is alone on the LHS, it is considered as the variable set by this equation
  log(dummy2) = k + 2;
  [ name = 'eq:dummy3', bar = 'baz' ]
  c(+1) = c;
end;

model_options(block);

model_remove('eq:dummy1', foo = 'eq:dummy2');

model_replace('ressource constraint', [ name = 'eq:dummy3', bar = 'baz' ]);
  c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
  c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
end;

var_remove dummy3;

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

S = load('../deterministic_simulations/ramst/Output/ramst_results.mat');
if any(size(oo_.endo_simul) ~= size(S.oo_.endo_simul)) || any(any(abs(oo_.endo_simul - S.oo_.endo_simul) > 1e-10))
  error('Model editing failure')
end
