// Test case with several exogenous variables (regression test for #1883)

var c k;
varexo x y;
parameters alph gam delt bet aa;

alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;

model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1) +y;
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
end;

initval;
y=0;
x = 1;
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k;
end;

steady;

check;

shocks(learnt_in=2);
  var x;
  periods 2:3;
  values 1.2;
end;

perfect_foresight_with_expectation_errors_setup(periods = 7);
perfect_foresight_with_expectation_errors_solver;

rplot c;
