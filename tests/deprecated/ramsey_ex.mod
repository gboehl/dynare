/* Tests the deprecated ramsey_policy command.
 *
 * The example is taken from Juillard, Michel (2011): User manual for optimal policy package, 
 * MONFISPOL FP7 project SSH-225149, Deliverable 1.1.2
*/

var pai, c, n, r, a;
varexo u;
parameters beta, rho, epsilon, omega, phi, gamma;

beta=0.99;
gamma=3;
omega=17;
epsilon=8;
phi=1;
rho=0.95;

model;
a = rho*a(-1)+u;
1/c = beta*r/(c(+1)*pai(+1));
pai*(pai-1)/c = beta*pai(+1)*(pai(+1)-1)/c(+1)+epsilon*phi*n^(gamma+1)/omega -exp(a)*n*(epsilon-1)/(omega*c);
exp(a)*n = c+(omega/2)*(pai-1)^2;
end;

initval;
r=1;
end;

initval;
a = 0;
pai = beta;
c = 0.9665;
n = 0.9673;
end;

shocks;
var u; stderr 0.008;
end;

planner_objective(ln(c)-phi*((n^(1+gamma))/(1+gamma)));
ramsey_policy(planner_discount=0.99,instruments=(r),order=1,nograph);

benchmark = load('../optimal_policy/Ramsey/ramsey_ex_initval/Output/ramsey_ex_initval_results.mat');

if any( [ max(abs(benchmark.oo_.steady_state-oo_.steady_state))>1e-5, ...
    max(abs(benchmark.oo_.dr.ys-oo_.dr.ys))>1e-5, ...
    max(max(abs(benchmark.oo_.dr.ghx-oo_.dr.ghx)))>1e-5, ...
    max(max(abs(benchmark.oo_.dr.ghu-oo_.dr.ghu)))>1e-5, ...
    max(max(abs(benchmark.oo_.dr.Gy-oo_.dr.Gy)))>1e-5, ...
    abs(benchmark.oo_.planner_objective_value.unconditional-oo_.planner_objective_value.unconditional)>1e-5, ...
    abs(benchmark.oo_.planner_objective_value.conditional.zero_initial_multiplier-oo_.planner_objective_value.conditional.zero_initial_multiplier)>1e-5, ...
    abs(benchmark.oo_.planner_objective_value.conditional.steady_initial_multiplier-oo_.planner_objective_value.conditional.steady_initial_multiplier)>1e-5] ) 
    error('ramsey_policy gives results inconsistent with ramsey_model+stoch_simul+evaluate_planner_objective')
end
