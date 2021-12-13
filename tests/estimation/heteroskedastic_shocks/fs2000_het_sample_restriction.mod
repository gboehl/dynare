// See fs2000.mod in the examples/ directory for details on the model
// tests heteroskedastic filter/smoother
// includes lagged exogenous variable introduced by preprocessor

@#include "fs2000_het_model.inc" 

shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
end;

steady;

check;
stoch_simul(order=1,loglinear);
forecast;

options_.solve_tolf = 1e-12;

heteroskedastic_shocks;
  var e_b;
  periods 100:120;
  values 0.01;

  var e_a;
  periods 100:120;
  scales 0;
end;

estimation(order=1,datafile='../fsdat_simul',first_obs=10,nobs=182,mode_compute=5,loglinear,mh_replic=0,smoother,filtered_vars,forecast=8,filter_step_ahead=[1:8],consider_all_endogenous,heteroskedastic_filter);

if M_.heteroskedastic_shocks.Qscale(strmatch('e_a',M_.exo_names,'exact'),91)~=0 && M_.heteroskedastic_shocks.Qscale(strmatch('e_b',M_.exo_names,'exact'),91)~=0.01
    error('first_obs is incorrectly handled.')
end