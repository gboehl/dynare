/* This model has a solve forward simple block that contains a recursive variable (a)
 and a feedback one (b).
 Regression test for bug fixed in preprocessor@740ea833f6b4a93c260b32e62f4302483af54f7a */

var c k a b;
varexo x e u;

parameters alph gam delt bet aa rho tau;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;
rho   = 0.95;
tau   = 0.025;


model(block, mfs=1);
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
a = tau*b + e;
b = tau*a+rho*b(-1) + u;
end;

initval;
x = 1;
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k;
end;

steady;

check;

shocks;
  var x;
  periods 1;
  values 1.2;
  var u;
  periods 1;
  values 1;
end;

model_info;

perfect_foresight_setup(periods=200);
perfect_foresight_solver;

rplot a;
