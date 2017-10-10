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
EfficiencyGrowth_ss = 1.02;
PopulationGrowth_ss = 1.02;

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

// Set initial condition for the forecast (we simulate a sample, forecast will start at the end of the sample)

histval;
    Efficiency(0) = 1;
    EfficiencyGrowth(0) = .5;
    Population(0) = 1;
    PopulationGrowth(0) = .5;
    PhysicalCapitalStock(0) = 1;
end;

oo__ = simul_backward_nonlinear_model([], 10, options_, M_, oo_);

initialcondition = dseries(transpose(oo__.endo_simul(:,10)), 2017Q1, M_.endo_names);


/* REMARKS
**
** dseries object initialcondition may have more variable than required.
*/

// Define the variables for which we want to compute the IRFs
listofvariables = {'Efficiency', 'Output'};

/* REMARKS
**
** If possible do not forecasts all the endogenous variables (this is the default).
*/

// Compute the forecasts
forecasts = backward_model_forecast(initialcondition, listofvariables, 16, true);

/* REMARKS
** 
** - The third argument is the number of periods for the forecast. Default is 8.
** - The last argument is for computing (or not) the confidence bands. If false (default), only point forecast is produced.
** - forecasts is a structure, each field is a dseries object for the point forecast (ie forecasts without innovations in the future),
** the mean forecast, the median forecast, the standard deviation of the predictive distribution and the lower/upper bounds of the 
** interval containing 95% of the predictive distribution.
*/

// Plot forecast for output
plot(forecasts.meanforecast.Output, '-k', 'linewidth', 2)
hold on
plot(forecasts.lb.Output,'-r')
plot(forecasts.ub.Output,'-r')
hold off

/* REMARKS
** 
** In this model there is no steady state (only a stable balanced growth paths), this explains the 
** shape of the forecast for Output.
*/