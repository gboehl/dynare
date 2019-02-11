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

histval;
  Efficiency(0) = 1;
  EfficiencyGrowth(0) = .5;
  Population(0) = 1;
  PopulationGrowth(0) = .5;
  PhysicalCapitalStock(0) = 1;
end;


/*
** First approach: Use the covariance matrix to define the impulses.
*/

// Define the shocks for which we want to compute the IRFs
listofshocks = {'e_x', 'e_n'};

// Define the variables for which we want to compute the IRFs
listofvariables = {'Efficiency', 'Population', 'Output'};

// Compute the IRFs
irfs1 = backward_model_irf([], dseries(), listofshocks, listofvariables, 50);

/*
** Second approach: Explicitely comunnicate the paths for the exogneous variables.
*/

x0 = zeros(50,2);
x0_1 = x0;
x0_2 = x0;
x0_1(1,1) = sqrt(M_.Sigma_e(1,1));
x0_2(1,2) = sqrt(M_.Sigma_e(2,2));

simshocks = {dseries(x0_1, 2, M_.exo_names), dseries(x0_2, 2, M_.exo_names) };

// Compute the IRFs
irfs2 = backward_model_irf([], dseries(), simshocks, listofvariables, 50);

/*
** Check that the two approaches do provide the same results.
*/

verbatim;
  if max(abs(irfs1.e_x.Efficiency.data-irfs2.experiment_1.Efficiency.data))>1e-12
    error('There is something wrong in the IRFs.')
  end
  if max(abs(irfs1.e_x.Output.data-irfs2.experiment_1.Output.data))>1e-12
    error('There is something wrong in the IRFs.')
  end
end;