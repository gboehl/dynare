// Test the [static]/[dynamic] tags, in the context of a block decomposed model

var c k;
varexo x;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;


model(block);
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
[dynamic] c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
[static] k = ((delt+bet)/(x*aa*alph))^(1/(alph-1));
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
end;

perfect_foresight_setup(periods=200);
perfect_foresight_solver;

rplot c;
rplot k;
