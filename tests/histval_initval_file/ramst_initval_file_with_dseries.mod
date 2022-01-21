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
(c/c(-2))^(-gam) - (1+bet)^(-1)*(aa*alph*x(+2)*k^(alph-1) + 1 - delt)*(c(+1)/c(-1))^(-gam);
end;

kstar = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
cstar = aa*kstar^alph-delt*kstar;

initval;
x = 1;
k = kstar;
c = cstar;
end;

steady;

// Instantiate dseries object with paths for the endogenous and exogenous variables
baseline= dseries([.81*kstar, .81*cstar, .81; .9*kstar, .9*cstar, .9; repmat([kstar, cstar, 1], 998, 1)], 2000Q1, {'k', 'c', 'x'});

initval_file(series=baseline, first_simulation_period=2000Q3, last_simulation_period=2050Q2);

perfect_foresight_setup;
perfect_foresight_solver;

Simulated_time_series(2000Q3:2050Q2)
