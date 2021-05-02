// See fs2000.mod in the examples/ directory for details on the model
// Tests that setting scale and value in same period for same shock is correctly filtered out

@#include "fs2000_het_model.inc" 

shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
end;

steady;

check;

options_.solve_tolf = 1e-12;

heteroskedastic_shocks;
  var e_b;
  periods 100:120;
  values 0.01;

  var e_a;
  periods 100:120;
  scales 0;

  // Wrongly set scale on top of value for the same shock and period
  var e_a;
  periods 110;
  values 0.01;
end;

estimation(order=1,datafile='../fsdat_simul',nobs=192,loglinear,mh_replic=0,smoother,filtered_vars,forecast=8,filter_step_ahead=[1:8],consider_all_endogenous,heteroskedastic_filter);
