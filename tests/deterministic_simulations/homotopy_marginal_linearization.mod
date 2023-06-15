// Example that triggers homotopy in perfect foresight simulation.
// Tests the homotopy_marginal_linearization_fallback and homotopy_step_size_increase_success_count options

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

perfect_foresight_solver(homotopy_max_completion_share = 0.7,
                         homotopy_step_size_increase_success_count = 3,
                         homotopy_marginal_linearization_fallback,
                         steady_solve_algo = 13);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end
