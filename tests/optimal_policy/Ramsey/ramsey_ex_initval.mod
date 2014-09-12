/* Mod file tests the functioning of the Ramsey command when used together with an initval-block
 * - Tests whether subsequent calls to ramsey_policy and stoch_simul are possible
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
ramsey_policy(planner_discount=0.99,order=1,instruments=(r));
oo_ramsey_policy_initval=oo_;
%save oo_ramsey_policy_initval oo_ramsey_policy_initval
stoch_simul(periods=0, order=1, irf=25, nograph);
if max(abs((oo_ramsey_policy_initval.steady_state-oo_.steady_state)))>1e-5 ...
    || max(abs(oo_ramsey_policy_initval.dr.ys-oo_.dr.ys))>1e-5 ...
    || max(max(abs(oo_ramsey_policy_initval.dr.ghx-oo_.dr.ghx)))>1e-5 ...
    || max(max(abs(oo_ramsey_policy_initval.dr.ghu-oo_.dr.ghu)))>1e-5 ...
    || max(max(abs(oo_ramsey_policy_initval.dr.Gy-oo_.dr.Gy)))>1e-5 ...
    || max(abs((oo_ramsey_policy_initval.planner_objective_value-oo_.planner_objective_value)))>1e-5
    error('Running stoch_simul after ramsey_policy leads to inconsistent results')
end
