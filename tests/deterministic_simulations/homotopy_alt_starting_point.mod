// Test for the homotopy_alt_starting_point option of perfect_foresight_solver

var Consumption, Capital, LoggedProductivity;

varexo LoggedProductivityInnovation;

parameters beta, alpha, delta, rho;

beta = .985;
alpha = 1/3;
delta = alpha/10;
rho = .9;

model;

[name='Euler equation']
1/Consumption = beta/Consumption(1)*(alpha*exp(LoggedProductivity(1))*Capital^(alpha-1)+1-delta);

[name='Physical capital stock law of motion']
Capital = exp(LoggedProductivity)*Capital(-1)^alpha+(1-delta)*Capital(-1)-Consumption;

[name='Logged productivity law of motion']
LoggedProductivity = rho*LoggedProductivity(-1)+LoggedProductivityInnovation;

end;

steady_state_model;
  LoggedProductivity = LoggedProductivityInnovation/(1-rho);
  Capital = (exp(LoggedProductivity)*alpha/(1/beta-1+delta))^(1/(1-alpha));
  Consumption = exp(LoggedProductivity)*Capital^alpha-delta*Capital;
end;

initval;
  LoggedProductivityInnovation = 0;
end;

steady;

endval;
  LoggedProductivityInnovation = 0.4;
end;

steady;

perfect_foresight_setup(periods=200);
perfect_foresight_solver(homotopy_alt_starting_point);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end