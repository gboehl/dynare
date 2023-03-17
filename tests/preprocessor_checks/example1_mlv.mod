// Tests for model local variables
// (including in params derivs file, i.e. with identification, see Dynare/preprocessor#13)

var c, y, k, a, h, b;
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

/* The following statement is a regression test for #1782.
   Here the “foo” variable definition depends on “bar”, but the symbol ID of
   “foo” will be smaller than the symbol ID of “bar”. */
model_local_variable foo $\text{foo}$;

model;
#bar = exp(b)*c;
#foo = bar/(exp(b(+1))*c(+1));
c*theta*h^(1+psi)=(1-alpha)*y;
k = beta*(foo
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

steady;
check;
stoch_simul;

varobs c;

identification(parameter_set=calibration);
