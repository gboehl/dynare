// Base simulation file adapted from ramst.mod

var c k;
varexo x;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;

model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1); // Resource constraint
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam); // Euler equation
end;

initval;
c = 1.2;
k = 12;
x = 1;
end;

shocks;
var x;
periods 2;
values 0.9;
end;

simul(periods=200,maxit=100);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end