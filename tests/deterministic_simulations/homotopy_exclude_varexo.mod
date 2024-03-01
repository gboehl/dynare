// Example that triggers homotopy in perfect foresight simulation.
// Tests the homotopy_exclude_varexo option

var Consumption, Capital, LoggedProductivity;

varexo LoggedProductivityInnovation dummy;

parameters beta, alpha, delta, rho;

beta = .985;
alpha = 1/3;
delta = alpha/10;
rho = .9;

model;
  1/Consumption = beta/Consumption(1)*(alpha*exp(LoggedProductivity(1))*Capital^(alpha-1)+1-delta)*dummy;
  Capital = exp(LoggedProductivity)*Capital(-1)^alpha+(1-delta)*Capital(-1)-Consumption;
  LoggedProductivity = rho*LoggedProductivity(-1)+LoggedProductivityInnovation;
end;

initval;
  LoggedProductivityInnovation = 0;
  dummy = 1;
end;

steady;

endval;
  LoggedProductivityInnovation = 1;
  dummy = 1;
  Consumption = 0.1;
  Capital = 1;
end;

perfect_foresight_setup(periods=200, endval_steady);

perfect_foresight_solver(homotopy_exclude_varexo = (dummy),
                         homotopy_max_completion_share = 0.7,
                         homotopy_linearization_fallback,
                         steady_solve_algo = 13);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

if oo_.deterministic_simulation.sim1.exo_simul(end, 1) ~= 0.7 ... % productivity must be rescaled
   || oo_.deterministic_simulation.sim1.exo_simul(end, 2) ~= 1    % but not dummy
    error('Option homotopy_exclude_varexo not working properly')
end
