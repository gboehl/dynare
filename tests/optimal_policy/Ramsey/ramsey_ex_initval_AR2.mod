/* Test the correctness of the Ramsey command when a Lagrange multiplier
 * appears with a lead ⩾2, and thus ends up in the definition of an auxiliary variable.
 *
 * This is related to issues #633, #1119 and #1133
 *
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
ramsey_model(planner_discount=0.99,instruments=(r));
steady;
stoch_simul(order=1, nograph);
evaluate_planner_objective;
