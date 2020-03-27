/*
  Tests that local_state_space_iteration_2 and local_state_space_iteration_k
  (for k=2) return the same results.
*/

var y, k, a, h, b, c;
varexo e, u;

parameters beta, rho, alpha, delta, theta, psi, tau;

alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;

phi   = 0.1;

model;
c*theta*h^(1+psi)=(1-alpha)*y;
k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
    *(exp(b(+1))*alpha*y(+1)+(1-delta)*k));
y = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
k = exp(b)*(y-c)+(1-delta)*k(-1);
a = rho*a(-1)+tau*b(-1) + e;
b = tau*a(-1)+rho*b(-1) + u;
end;

initval;
y = 1.08068253095672;
c = 0.80359242014163;
h = 0.29175631001732;
k = 11.08360443260358;
a = 0;
b = 0;
e = 0;
u = 0;
end;

shocks;
var e; stderr 0.009;
var u; stderr 0.009;
var e, u = phi*0.009*0.009;
end;

stoch_simul(order=2, irf=0, k_order_solver);

nparticles = 100;

/* We generate particles using realistic distributions (though this is not
   strictly needed) */
state_idx = oo_.dr.order_var((M_.nstatic+1):(M_.nstatic+M_.npred+M_.nboth));
yhat = chol(oo_.var(state_idx,state_idx))*randn(M_.npred+M_.nboth, nparticles);
epsilon = chol(M_.Sigma_e)*randn(M_.exo_nbr, nparticles);

dr = oo_.dr;

tStart1 = tic; for i=1:10000, ynext1 = local_state_space_iteration_2(yhat, epsilon, dr.ghx, dr.ghu, dr.ys(dr.order_var)+0.5*dr.ghs2, dr.ghxx, dr.ghuu, dr.ghxu, 1); end, tElapsed1 = toc(tStart1);

tStart2 = tic; for i=1:10000, ynext2 = local_state_space_iteration_k(yhat, epsilon, dr, M_, options_); end, tElapsed2 = toc(tStart2);

if max(max(abs(ynext1-ynext2))) > 1e-14
    error('Inconsistency between local_state_space_iteration_2 and local_state_space_iteration_k')
end

if tElapsed1<tElapsed2
    skipline()
    dprintf('local_state_space_iteration_2 is %5.2f times faster than local_state_space_iteration_k', tElapsed2/tElapsed1)
    skipline()
else
    skipline()
    dprintf('local_state_space_iteration_2 is %5.2f times slower than local_state_space_iteration_k', tElapsed1/tElapsed2)
    skipline()
end
