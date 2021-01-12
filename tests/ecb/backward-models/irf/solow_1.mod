close all

var Efficiency $A$
    EfficiencyGrowth $X$
    Population $L$
    PopulationGrowth $N$
    Output $Y$
    PhysicalCapitalStock $K$ ;

varexo e_x $\varepsilon_x$
       e_n $\varepsilon_n$;

parameters alpha $\alpha$
	   delta $\delta$
	   s $s$
           rho_x $\rho_x$
           rho_n $\rho_n$
           EfficiencyGrowth_ss $X^{\star}$
           PopulationGrowth_ss $N^{\star}$ ;

alpha = .33;
delta = .02;
s     = .20;
rho_x = .90;
rho_n = .95;
EfficiencyGrowth_ss = 1.0;
PopulationGrowth_ss = 1.0;

model;
    Efficiency = EfficiencyGrowth*Efficiency(-1);
    EfficiencyGrowth/EfficiencyGrowth_ss = (EfficiencyGrowth(-1)/EfficiencyGrowth_ss)^(rho_x)*exp(e_x);
    Population = PopulationGrowth*Population(-1);
    PopulationGrowth/PopulationGrowth_ss = (PopulationGrowth(-1)/PopulationGrowth_ss)^(rho_n)*exp(e_n);
    Output = PhysicalCapitalStock(-1)^alpha*(Efficiency*Population)^(1-alpha);
    PhysicalCapitalStock = (1-delta)*PhysicalCapitalStock(-1) + s*Output;
end;


shocks;
    var e_x = 0.005;
    var e_n = 0.001;
end;

/*
** FIRST APPROACH (PASS A DSERIES OBJECT FOR THE INITIAL CONDITION)
*/

// Define a dseries object.
ds = dseries(repmat([1, .5, 1, .5, 0, 1], 10, 1), 1990Q1, M_.endo_names, M_.endo_names_tex);

// Alternatively we could build this object with a stochastic simulation of the model.
/*
ds = simul_backward_model(dseries([1, .5, 1, .5, 0, 1], 1990Q1, cellstr(M_.endo_names)), 10);
names=regexp(ds.name, 'e_\w*');
idxs = [];
for j=1:length(names)
    if isempty(names{j})
        idxs = [idxs j];
    end
end
ds = ds{idxs};
*/

// Set the initial condition for the IRFs using observation 1991Q1 in ds.
set_historical_values ds 1991Q1;

// Define the shocks for which we want to compute the IRFs
listofshocks = {'e_x', 'e_n'};

// Define the variables for which we want to compute the IRFs
listofvariables = {'Efficiency', 'Population', 'Output'};

// Compute the IRFs
irfs = backward_model_irf(ds, dseries(), listofshocks, listofvariables, 50);  // 10 is the number of periods (default value is 40).

// Plot an IRF (shock on technology)
figure(1)
plot(irfs.e_x.Output)
legend('Output')
title('IRF (shock on the growth rate of efficiency)')

/*
** SECOND APPROACH (USE THE HISTVAL BLOCK FOR SETTING THE INITIAL CONDITION)
**
** In this case the first argument of the backward_model_irf routine is a date object (for specifying
** the period of the initial condition).
*/

/*
histval;
  Efficiency(0) = 1;
  EfficiencyGrowth(0) = .5;
  Population(0) = 1;
  PopulationGrowth(0) = .5;
  PhysicalCapitalStock(0) = 1;
end;
*/
// Define the shocks for which we want to compute the IRFs
listofshocks = {'e_x', 'e_n'};

// Define the variables for which we want to compute the IRFs
listofvariables = {'Efficiency', 'Population', 'Output'};

// Compute the IRFs
irfs = backward_model_irf([], dseries(), listofshocks, listofvariables, 50);  // 10 is the number of periods (default value is 40).

// Plot an IRF (shock on technology)
figure(2)
plot(irfs.e_x.Output)
legend('Output')
title('IRF (shock on the growth rate of efficiency)')