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
var e = .01;
var u = .02;
end;

stoch_simul(periods=200);

epilogue;
g = log((y/y(-4))^(1/4));
x = .7*x(-1)+g(-1)*alpha;
end;

ds = dseries(oo_.endo_simul', 2000Q1, M_.endo_names);
ds = [ds dseries(randn(7,1), 2000Q1, 'x')];
ds = example1.epilogue(M_.params, ds);
