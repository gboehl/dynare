% Test whether preprocessor recognizes state variables introduced by optimal policy Github #1193
        
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

histval;
r(0)=1;
end;

steady_state_model;
a = 0;
pai = beta*r;
c = find_c(0.96,pai,beta,epsilon,phi,gamma,omega);
n = c+(omega/2)*(pai-1)^2;
end;

shocks;
var u; stderr 0.008;
var u;
periods 1;
values 1;
end;
options_.dr_display_tol=0;
planner_objective(ln(c)-phi*((n^(1+gamma))/(1+gamma)));
ramsey_policy(planner_discount=0.99,order=1,instruments=(r),periods=500);