// Example 1 from Collard's guide to Dynare
var y, c, k, a, h, b;
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

external_function(nargs=2, name=extFunNoDerivs);

model(use_dll);
c*theta*h^(1+psi)=(1-alpha)*y;
k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
    *(exp(b(+1))*alpha*y(+1)+(1-delta)*k));
y = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
k = exp(b)*(y-c)+extFunNoDerivs((1-delta),k(-1));
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

stoch_simul;

L = load(['benchmark' filesep 'Output' filesep 'benchmark_results.mat']);
if max(max(abs(L.oo_.dr.ghu - oo_.dr.ghu))) > 1e-9
  error('Failure in external function')
end
if max(max(abs(L.oo_.dr.ghx - oo_.dr.ghx))) > 1e-9
  error('Failure in external function')
end
if max(max(abs(L.oo_.dr.ghxu - oo_.dr.ghxu))) > 1e-4
  error('Failure in external function')
end
if max(max(abs(L.oo_.dr.ghxx - oo_.dr.ghxx))) > 2e-3
  error('Failure in external function')
end
if max(max(abs(L.oo_.dr.ghuu - oo_.dr.ghuu))) > 1e-5
  error('Failure in external function')
end
if max(max(abs(L.oo_.dr.ghs2 - oo_.dr.ghs2))) > 1e-9
  error('Failure in external function')
end
if max(max(abs(L.oo_.var - oo_.var))) > 1e-12
  error('Failure in external function')
end
