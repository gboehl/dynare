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

@#define fake_backward = 1

var y1 y2 y3 ;

@#if fake_backward
  var z ;
@#endif

varexo e1 e2 e3 ;

parameters a b c ;

a = .1;
b = .2;
c = .3;

model;
  @#if fake_backward
    z = a*z(-1) ;
  @#endif
  y1 = y2*y3+exp(e1) ;
  y2 = 1 + y1*exp(e2) ;
  y3 = e3 ;
end;

@#if fake_backward
  histval;
    z(0) = .0 ;
  end;
@#endif

shocks;
    var e1 = .01;
    var e2 = .01;
    var e3 = .02;
end;

TrueData = simul_backward_model([], 1000);

// Set the periods where some of the endogenous variables will be constrained.
subsample = 2Y:101Y;

// Load the generated data
SimulatedData = copy(TrueData);

// Set the constrained paths for the endogenous variables.
constrainedpaths = SimulatedData{'y1','y2','y3'}(subsample);

// Set the instruments (innovations used to control the paths for the endogenous variables).
exogenousvariables = dseries([NaN(100, 3)], '2Y', M_.exo_names);

/* REMARK
**
** Here we will control y1, y2, and y3 with  e1, e2 and e3.
**
*/

// Invert the model by calling the model_inversion routine.
options_.dynatol.f = 1e-9;
[endogenousvariables, exogenousvariables] = model_inversion(constrainedpaths, exogenousvariables, SimulatedData, M_, options_, oo_);

// Check that all the constraints are satisfied.
if max(abs(constrainedpaths(subsample).y1.data-endogenousvariables(subsample).y1.data))>1e-12
   error('Constraint on y1 path is not satisfied!')
end

if max(abs(constrainedpaths(subsample).y2.data-endogenousvariables(subsample).y2.data))>1e-12
   error('Constraint on y2 path is not satisfied!')
end

if max(abs(constrainedpaths(subsample).y3.data-endogenousvariables(subsample).y3.data))>1e-12
   error('Constraint on y3 path is not satisfied!')
end

// Check consistency of the results.
if max(abs(exogenousvariables(subsample).e1.data-SimulatedData(subsample).e1.data))>1e-8
   error('Model inversion is not consistent with true innovations (e1)')
end

if max(abs(exogenousvariables(subsample).e2.data-SimulatedData(subsample).e2.data))>1e-8
   error('Model inversion is not consistent with true innovations (e2)')
end

if max(abs(exogenousvariables(subsample).e3.data-SimulatedData(subsample).e3.data))>1e-8
   error('Model inversion is not consistent with true innovations (e3)')
end
