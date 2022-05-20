// --+ options: stochastic +--

/* © 2022 Dynare Team
 *
 * This file is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * It is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with the file.  If not, see <http://www.gnu.org/licenses/>.
 */

var Efficiency                  // $A$
    EfficiencyGrowth            // $X$
    Population                  // $L$
    PopulationGrowth            // $N$
    Output                      // $Y$
    PhysicalCapitalStock ;      // $K$

varexo e_x   // $\varepsilon_x$
       e_n ; // $\varepsilon_n$

parameters alpha                               // $\alpha$
           epsilon                             // $\varepsilon$
	   delta                               // $\delta$
	   s                                   // $s$
           rho_x                               // $\rho_x$
           rho_n                               // $\rho_n$
           EfficiencyGrowth_ss                 // $X^{\star}$
           PopulationGrowth_ss ;               // $N^{\star}$

alpha   = .33;
epsilon = .70;
delta   = .02;
s       = .20;
rho_x   = .90;
rho_n   = .95;
EfficiencyGrowth_ss = 1.02;
PopulationGrowth_ss = 1.02;

if s>delta*alpha^(-epsilon/(epsilon-1))
    disp('The model admits a unique positive steady state.')
end


model;
    Efficiency = EfficiencyGrowth*Efficiency(-1);
    EfficiencyGrowth/EfficiencyGrowth_ss = (EfficiencyGrowth(-1)/EfficiencyGrowth_ss)^(rho_x)*exp(e_x);
    Population = PopulationGrowth*Population(-1);
    PopulationGrowth/PopulationGrowth_ss = (PopulationGrowth(-1)/PopulationGrowth_ss)^(rho_n)*exp(e_n);
    Output = (alpha*PhysicalCapitalStock(-1)^((epsilon-1)/epsilon)+(1-alpha)*(Efficiency*Population)^((epsilon-1)/epsilon))^(epsilon/(epsilon-1));
    PhysicalCapitalStock = (1-delta)*PhysicalCapitalStock(-1) + s*Output;
end;

histval;
    Efficiency(0) = 1;
    EfficiencyGrowth(0) = 1.02;
    Population(0) = 1;
    PopulationGrowth(0) = 1.02;
    PhysicalCapitalStock(0) = 1;
end;

shocks;
    var e_x = 0.005;
    var e_n = 0.001;
end;

TrueData = simul_backward_nonlinear_model([], 100, options_, M_, oo_);

// Set the periods where some of the endogenous variables will be constrained.
subsample = 2Y:100Y;

// Copy the generated data
SimulatedData = copy(TrueData);

// Set the constrained paths for the endogenous variables (Output and PhysicalCapitalStock).
constrainedpaths = SimulatedData{'Output','PhysicalCapitalStock'}(subsample);

// Set the instruments (innovations used to control the paths for the endogenous variables).
exogenousvariables = dseries(NaN(100, 2), 1Y, {'e_x';'e_n'});

// Invert the model by calling the model_inversion routine.
[endogenousvariables, exogenousvariables] = model_inversion(constrainedpaths, exogenousvariables, SimulatedData, M_, options_, oo_);

// Check that all the constraints are satisfied.
if max(abs(constrainedpaths.Output(subsample).data-endogenousvariables.Output(subsample).data))>1e-12
   error('Constraint on Output path is not satisfied!')
end

if max(abs(constrainedpaths.PhysicalCapitalStock(subsample).data-endogenousvariables.PhysicalCapitalStock(subsample).data))>1e-12
   error('Constraint on PhysicalCapitalStock path is not satisfied!')
end

// Check the consistency of the deduced innovations with the true innovations.

/*  REMARK
**
**  In this model the two innovatons, on population growth and efficiency growth, cannot be identified simulatneously because
**  what matters for the dynamic of the physical capital stock and the output is the same non linear function of the two innovations.
**  This explains why we only check that the product of Efficiency and Population are the same in the simulated (true) data and in the
**  data returned by the inversion routine.
*/

TrueEfficiency = TrueData.Efficiency;
TruePopulation = TrueData.Population;

TrueEfficiencyTimesPopulation = TrueEfficiency*TruePopulation;
EfficiencyTimesPopulation = endogenousvariables.Efficiency*endogenousvariables.Population;

if max(abs(TrueEfficiencyTimesPopulation(2Y:100Y).data-EfficiencyTimesPopulation(2Y:100Y).data))>1e-5
   error('Model inversion is not consitent with true innovations.')
end
