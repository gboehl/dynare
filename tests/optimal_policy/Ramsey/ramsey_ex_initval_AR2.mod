/* Mod file tests the correctness of the Ramsey command when Auxiliary variables for AR2 are used
 * The example is adapted from Juillard, Michel (2011): User manual for optimal policy package, 
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
log(a) = rho*log(a(-2))+u;
1/c = beta*r/(c(+1)*pai(+1));
pai*(pai-1)/c = beta*pai(+1)*(pai(+1)-1)/c(+1)+epsilon*phi*n^(gamma+1)/omega -a*n*(epsilon-1)/(omega*c);
a*n = c+(omega/2)*(pai-1)^2;
end;

initval;
r=1;
end;

initval;
a = 1;
pai = beta;
c = 0.9665;
n = 0.9673;
end;

shocks;
var u; stderr 0.008;
end;

planner_objective(ln(c)-phi*((n^(1+gamma))/(1+gamma)));
ramsey_policy(planner_discount=0.99,order=1,instruments=(r));
oo_ramsey_policy_initval_AR2=oo_;
stoch_simul(periods=0, order=1, irf=25, nograph);
if max(abs((oo_ramsey_policy_initval_AR2.steady_state-oo_.steady_state)))>1e-5 ...
    || max(abs(oo_ramsey_policy_initval_AR2.dr.ys-oo_.dr.ys))>1e-5 ...
    || max(max(abs(oo_ramsey_policy_initval_AR2.dr.ghx-oo_.dr.ghx)))>1e-5 ...
    || max(max(abs(oo_ramsey_policy_initval_AR2.dr.ghu-oo_.dr.ghu)))>1e-5 ...
    || max(max(abs(oo_ramsey_policy_initval_AR2.dr.Gy-oo_.dr.Gy)))>1e-5 ...
    || max(abs((oo_ramsey_policy_initval_AR2.planner_objective_value-oo_.planner_objective_value)))>1e-5
    error('Running stoch_simul after ramsey_policy leads to inconsistent results')
end