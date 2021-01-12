// --+ options: json=compute +--

var y, c, k, h;
var(rows = 2) A; // A == [a; b]
varexo(rows = 2) E; // E == [e; u]
parameters beta, alpha, delta, theta, psi;
parameters(rows=2,cols=2) P;

alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;

phi   = 0.1;

P = [rho, tau; tau/2, rho/3];

external_function(nargs=3, name=extFunNoDerivsMatrix, in_arg_dim=[2,2,2],
                  out_arg_dim=2, first_deriv_provided=extFunFirstDerivMatrix);

model(use_dll);
c*theta*h^(1+psi)=(1-alpha)*y;
k = beta*(((exp(A[2])*c)/(exp(A[2](+1))*c(+1)))
          *(exp(A[2](+1))*alpha*y(+1)+(1-delta)*k));
y = exp(A[1])*(k(-1)^alpha)*(h^(1-alpha));
k = exp(A[2])*(y-c)+(1-delta)*k(-1) ;
A = extFunNoDerivsMatrix(P[:,1], transpose(P[:,2]), A(-1)) + E;
end;

initval;
y = 1.08068253095672;
c = 0.80359242014163;
h = 0.29175631001732;
k = 11.08360443260358;

A[1] = 0;
A[2] = 0;
E[:] = 0;

end;

shocks;
var E[:]; stderr 0.009;
var E[1], E[2] = phi*0.009*0.009;
end;

stoch_simul;

L = load('benchmark_results.mat');
if max(max(abs(L.oo_.dr.ghu - oo_.dr.ghu))) > 1e-12
  error('Failure in matrix notation with external function')
end
if max(max(abs(L.oo_.dr.ghx - oo_.dr.ghx))) > 1e-12
  error('Failure in matrix notation with external function')
end
if max(max(abs(L.oo_.dr.ghxu - oo_.dr.ghxu))) > 1e-12
  error('Failure in matrix notation with external function')
end
if max(max(abs(L.oo_.dr.ghxx - oo_.dr.ghxx))) > 1e-12
  error('Failure in matrix notation with external function')
end
if max(max(abs(L.oo_.dr.ghuu - oo_.dr.ghuu))) > 1e-12
  error('Failure in matrix notation with external function')
end
if max(max(abs(L.oo_.dr.ghs2 - oo_.dr.ghs2))) > 1e-12
  error('Failure in matrix notation with external function')
end
if max(max(abs(L.oo_.var - oo_.var))) > 1e-12
  error('Failure in matrix notation with external function')
end
