// Example 1 from Collard's guide to Dynare
var c, a, h, b;

verbatim;
% I want these comments included in
% example1.m 1999q1 1999y
%
var = 1;
end;

parameters beta, rho, theta, psi, tau;

rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
psi   = 0;
theta = 2.95;

phi   = 0.1;

model;
c*theta*h^(1+psi)=(1-alpha)*y;
k|e = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
    *(exp(b(+1))*alpha|p*y(+1)+(1-delta)*k));
y|e = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
k = exp(b)*(y-c)+(1-delta|p)*k(-1);
a = rho*a(-1)+tau*b(-1) + e|x;
b = tau*a(-1)+rho*b(-1) + u|x;
end;

alpha = 0.36;
delta = 0.025;

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
