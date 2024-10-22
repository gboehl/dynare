// check shocks on several periods
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
x = 1;
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k +1 ;
end;

steady;

check;

shocks;
var x;
periods 1 2 3 4;
values 1.1 1.2 1.3 1.4;
end;

perfect_foresight_setup(periods=200);
perfect_foresight_solver;

rplot c;
rplot k;
