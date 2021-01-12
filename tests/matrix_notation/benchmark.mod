/* This file is used as a benchmark against with the other tests with external
   functions are compared.
   It is almost the same as example1.mod, except that the shocks process is
   nonlinear (in order to have a non-zero Hessian of the
   external function) */

var y, c, k, h, a, b;
varexo e, u;
parameters beta, alpha, delta, theta, psi, rho, tau;

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
k = exp(b)*(y-c)+(1-delta)*k(-1) ;
a = rho*a(-1)^2+5*tau*b(-1)^2 + e;
b = 3*tau/2*a(-1)^2+7*rho/3*b(-1)^2 + u;
end;

initval;
y = 1.08068253095672;
c = 0.80359242014163;
h = 0.29175631001732;
k = 11.08360443260358;
end;

shocks;
var e; stderr 0.009;
var u; stderr 0.009;
var e, u = phi*0.009*0.009;
end;

stoch_simul;
