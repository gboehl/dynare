var Efficiency                  // $A$
    EfficiencyGrowth            // $X$
    Population                  // $L$
    PopulationGrowth            // $N$
    Output                      // $Y$
    PhysicalCapitalStock ;      // $K$

varexo e_x   // $\varepsilon_x$
       e_n ; // $\varepsilon_n$

parameters alpha                               // $\alpha$
	   delta                               // $\delta$
	   s                                   // $s$
           rho_x                               // $\rho_x$
           rho_n                               // $\rho_n$
           EfficiencyGrowth_ss                 // $X^{\star}$
           PopulationGrowth_ss ;               // $N^{\star}$

alpha = .33;
delta = .02;
s     = .20;
rho_x = .90;
rho_n = .95;
EfficiencyGrowth_ss = 1.00; // Do not change this calibration
PopulationGrowth_ss = 1.00; // Do not change this calibration

model;
    Efficiency = EfficiencyGrowth*Efficiency(-1);
    EfficiencyGrowth/EfficiencyGrowth_ss = (EfficiencyGrowth(-1)/EfficiencyGrowth_ss)^(rho_x)*exp(e_x);
    Population = PopulationGrowth*Population(-1);
    PopulationGrowth/PopulationGrowth_ss = (PopulationGrowth(-1)/PopulationGrowth_ss)^(rho_n)*exp(e_n);
    Output = PhysicalCapitalStock(-1)^alpha*(Efficiency*Population)^(1-alpha);
    PhysicalCapitalStock = (1-delta)*PhysicalCapitalStock(-1) + s*Output;
end;

d = dseries([.5 1 1 1.04 15], 2000Q1, {'Efficiency'; 'EfficiencyGrowth'; 'Population'; 'PopulationGrowth'; 'PhysicalCapitalStock'});

histval;
    Efficiency(0) = .5;
    EfficiencyGrowth(0) = 1;
    Population(0) = 1;
    PopulationGrowth(0) = 1.04;
    PhysicalCapitalStock(0) = 15;
end;

LongRunEfficiency = d.Efficiency*d.EfficiencyGrowth^(rho_x/(1-rho_x));
LongRunPopulation = d.Population*d.PopulationGrowth^(rho_n/(1-rho_n));
LongRunEfficiencyGrowth = EfficiencyGrowth_ss; 
LongRunPopulationGrowth = PopulationGrowth_ss;
LongRunIntensiveCapitalStock = (s/(LongRunEfficiencyGrowth*LongRunPopulationGrowth-1+delta))^(1/(1-alpha)); //LongRunEfficiencyGrowth*LongRunPopulationGrowth*

precision = 1e-6;
T = 5*floor(log(precision)/log(max(rho_x, rho_n)));

e = dseries(zeros(T, 2), 2000Q2, {'e_x'; 'e_n'});

simulations = simul_backward_model([], T, e);

if abs(simulations.Efficiency.data(end)-LongRunEfficiency.data)>1e-10
   error('Wrong long run level (Efficiency, diff=%s)!', num2str(abs(simulations.Efficiency.data(end)-LongRunEfficiency.data)))
end

if abs(simulations.Population.data(end)-LongRunPopulation.data)>1e-10
   error('Wrong long run level (Population, diff=%s)!', num2str(abs(simulations.Population.data(end)-LongRunPopulation.data)))
end

IntensiveCapitalStock = simulations.PhysicalCapitalStock/(simulations.Efficiency*simulations.Population);

if abs(IntensiveCapitalStock.data(end)-LongRunIntensiveCapitalStock)>1e-3
    error('Wrong long run level (Intensive capital stock)!')
end
