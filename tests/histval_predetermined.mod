/* Check that histval and predetermined_variables interact correctly, in the
   case where no explicit lag appears in the model block (though lags appear
   implicitly via the “predetermined_variables” statement).

   This is a regression test for preprocessor#47.
*/

var y, c, k, a, h, b;
varexo e, u;

predetermined_variables a b k;

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
k(+1) = beta*(((exp(b(+1))*c)/(exp(b(+2))*c(+1)))
        *(exp(b(+2))*alpha*y(+1)+(1-delta)*k(+1)));
y = exp(a(+1))*(k^alpha)*(h^(1-alpha));
k = exp(b(+1))*(y-c)+(1-delta)*k;
a(+1) = rho*a+ e;
b(+1) = rho*b+ u;
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

histval;
b(0) = 0.1;
a(0) = 0.3;
end;

stoch_simul(nograph, periods = 200);

forecast;

