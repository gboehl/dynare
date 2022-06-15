var Capital, Output, Labour, Consumption, Efficiency, efficiency, ExpectedTerm;

varexo EfficiencyInnovation;

parameters beta, theta, tau, alpha, psi, delta, rho, effstar, sigma2;

beta    =   0.9900;
theta   =   0.3570;
tau     =   2.0000;
alpha   =   0.4500;
psi     =  -0.1000;
delta   =   0.0200;
rho     =   0.8000;
effstar =   1.0000;
sigma2  =   0;

model;

  // Eq. n°1:
  efficiency = rho*efficiency(-1) + EfficiencyInnovation;

  // Eq. n°2:
  Efficiency = effstar*exp(efficiency);

  // Eq. n°3:
  Output = Efficiency*(alpha*(Capital(-1)^psi)+(1-alpha)*(Labour^psi))^(1/psi);

  // Eq. n°4:
  Capital = Output-Consumption + (1-delta)*Capital(-1);

  // Eq. n°5:
  ((1-theta)/theta)*(Consumption/(1-Labour)) - (1-alpha)*(Output/Labour)^(1-psi);

  // Eq. n°6:
  (((Consumption^theta)*((1-Labour)^(1-theta)))^(1-tau))/Consumption  = ExpectedTerm(1);

  // Eq. n°7:
  ExpectedTerm = beta*((((Consumption^theta)*((1-Labour)^(1-theta)))^(1-tau))/Consumption)*(alpha*((Output/Capital(-1))^(1-psi))+(1-delta));

end;

steady_state_model;
efficiency = EfficiencyInnovation/(1-rho);
Efficiency = effstar*exp(efficiency);
Output_per_unit_of_Capital=((1/beta-1+delta)/alpha)^(1/(1-psi));
Consumption_per_unit_of_Capital=Output_per_unit_of_Capital-delta;
Labour_per_unit_of_Capital=(((Output_per_unit_of_Capital/Efficiency)^psi-alpha)/(1-alpha))^(1/psi);
Output_per_unit_of_Labour=Output_per_unit_of_Capital/Labour_per_unit_of_Capital;
Consumption_per_unit_of_Labour=Consumption_per_unit_of_Capital/Labour_per_unit_of_Capital;

% Compute steady state share of capital.
ShareOfCapital=alpha/(alpha+(1-alpha)*Labour_per_unit_of_Capital^psi);

% Compute steady state of the endogenous variables.
Labour=1/(1+Consumption_per_unit_of_Labour/((1-alpha)*theta/(1-theta)*Output_per_unit_of_Labour^(1-psi)));
Consumption=Consumption_per_unit_of_Labour*Labour;
Capital=Labour/Labour_per_unit_of_Capital;
Output=Output_per_unit_of_Capital*Capital;
ExpectedTerm=beta*((((Consumption^theta)*((1-Labour)^(1-theta)))^(1-tau))/Consumption)
             *(alpha*((Output/Capital)^(1-psi))+1-delta);
end;

steady;

ik = varlist_indices('Capital',M_.endo_names);
CapitalSS = oo_.steady_state(ik);

histval;
Capital(0) = CapitalSS/2;
end;


perfect_foresight_setup(periods=400);
perfect_foresight_solver(stack_solve_algo=0);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

oo0 = oo_;

perfect_foresight_setup(periods=400);
perfect_foresight_solver(stack_solve_algo=1);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

oo1 = oo_;

maxabsdiff = max(max(abs(oo0.endo_simul-oo1.endo_simul)));

if max(max(abs(oo0.endo_simul-oo1.endo_simul)))>options_.dynatol.x
    error('stack_solve_algo={0,1} return different paths for the endogenous variables!')
else
    skipline()
    fprintf('Maximum (absolute) differrence between paths is %s', num2str(maxabsdiff))
    skipline()
end

% Also test stack_solve_algo=6, which is a synonymous for stack_solve_algo=1

perfect_foresight_setup(periods=400);
perfect_foresight_solver(stack_solve_algo=6);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

oo6 = oo_;

maxabsdiff = max(max(abs(oo0.endo_simul-oo6.endo_simul)));

if max(max(abs(oo0.endo_simul-oo6.endo_simul)))>options_.dynatol.x
    error('stack_solve_algo={0,6} return different paths for the endogenous variables!')
else
    skipline()
    fprintf('Maximum (absolute) differrence between paths is %s', num2str(maxabsdiff))
    skipline()
end
