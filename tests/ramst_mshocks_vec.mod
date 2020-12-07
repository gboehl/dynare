// Test mshocks block

var c k;
varexo x;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;


model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
end;

initval;
x = 2;
k = ((delt+bet)/(1.0*aa*x*alph))^(1/(alph-1));
c = aa*x*k^alph-delt*k;
end;

steady;

check;

xvalues = [1.2; 0.8; 0.8];

mshocks;
var x;
periods 1:3;
values (xvalues);
end;

perfect_foresight_setup(periods=200);
perfect_foresight_solver;

rplot c;
rplot k;

if ~all(oo_.exo_simul(M_.maximum_lag+(1:3)) == [ 2.4; 1.6; 1.6]) ...
    || ~all(oo_.exo_simul(M_.maximum_lag+(4:200)) == 2)
  error("mshocks not correctly applied")
end
