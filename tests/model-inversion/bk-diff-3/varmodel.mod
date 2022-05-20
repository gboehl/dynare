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

@#define USE_HISTVAL = 0

var y ;

varexo e ;

parameters rho ;

rho = 1.0;

model(linear);
    diff(y) =  rho*diff(y(-1)) + e ;
end;

shocks;
    var e = 1.0;
end;

steady;

check;

innovations = dseries(ones(200,1), 3Y, 'e');

@#if USE_HISTVAL
    histval;
        y(-1) =  0.0;
    	y(0) =  1.0;
    end;
    TrueData = simul_backward_model([], 200, innovations);
@#else
    initdata = zeros(2, 2);
    initdata(1,1) = 0.0;
    initdata(2,1) = 1.0;
    initialconditions = dseries(initdata, 1Y, vertcat(M_.endo_names(1:M_.orig_endo_nbr), M_.exo_names));
    TrueData = simul_backward_model(initialconditions, 200, innovations);
@#endif

y0 = TrueData(2Y).y.data;
yts = TrueData(3Y:10Y).y.data;
yex = 1+cumsum(y0+cumsum(ones(8,1)));

assert(~any(yts-yex), 'Simulation is wrong!')

// Set the periods where some of the endogenous variables will be constrained.
subsample = 3Y:100Y;

// Load the generated data
SimulatedData = copy(TrueData);

// Set the constrained paths for the endogenous variables (Output and PhysicalCapitalStock).
constrainedpaths = SimulatedData{'y'}(subsample);

// Set the instruments (innovations used to control the paths for the endogenous variables).
exogenousvariables = dseries(NaN(98, 1), 3Y, M_.exo_names);

// Invert the model by calling the model_inversion routine.
[endogenousvariables, exogenousvariables] = model_inversion(constrainedpaths, exogenousvariables, SimulatedData, M_, options_, oo_);

assert(all(exogenousvariables.e.data==1), 'Inversion is wrong!');
