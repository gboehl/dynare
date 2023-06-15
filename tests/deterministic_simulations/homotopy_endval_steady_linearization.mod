// Example that triggers homotopy in perfect foresight simulation.
// Tests the endval_steady, homotopy_linearization_fallback options, and a few more.

var Consumption, Capital, LoggedProductivity;

varexo LoggedProductivityInnovation;

parameters beta, alpha, delta, rho;

beta = .985;
alpha = 1/3;
delta = alpha/10;
rho = .9;

model;
  1/Consumption = beta/Consumption(1)*(alpha*exp(LoggedProductivity(1))*Capital^(alpha-1)+1-delta);
  Capital = exp(LoggedProductivity)*Capital(-1)^alpha+(1-delta)*Capital(-1)-Consumption;
  LoggedProductivity = rho*LoggedProductivity(-1)+LoggedProductivityInnovation;
end;

initval;
  LoggedProductivityInnovation = 0;
end;

steady;

endval;
  LoggedProductivityInnovation = 1;
  Consumption = 0.1;
  Capital = 1;
end;

perfect_foresight_setup(periods=200, endval_steady);

perfect_foresight_solver(homotopy_initial_step_size = 0.5,
                         homotopy_max_completion_share = 0.6,
                         homotopy_min_step_size = 0.0001,
                         homotopy_linearization_fallback,
                         steady_solve_algo = 13, steady_tolf=1e-7);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end
