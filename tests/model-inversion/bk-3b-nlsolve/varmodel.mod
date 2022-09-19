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

var y1 y2 y3 ;

varexo e1 e2 e3 ;

parameters a11 a12 a13 a21 a22 a23 a31 a32 a33 b11 b12 b13 b22 b23 b33 ;

/*
** Simulate the elements of the first order autoregressive matrix (we impose stability of the model, note that
** inversion fails if the model is explosive, ie the autoregressive matrix has at least one root greater than
** one in modulus) probably because of the propagation of roundoff errors.
*/

D = diag([.9 .7 .8]);
P = randn(3,3);
A = P*D*inv(P);

a11 = A(1,1);
a12 = A(1,2);
a13 = A(1,3);
a21 = A(2,1);
a22 = A(2,2);
a23 = A(2,3);
a31 = A(3,1);
a32 = A(3,2);
a33 = A(3,3);
b11 =  .10;
b12 = -.30;
b13 =  .05;
b22 =  .20;
b23 = -.05;
b33 =  .10;

model;
    y1 = a11*y1(-1) + a12*y2(-1) + a13*y3(-1) + b11*e1 + b12*e2 + b13*e3 ;
    y2 = a21*y1(-1) + a22*y2(-1) + a23*y3(-1)          + b22*e2 + b23*e3 ;
    y3 = a31*y1(-1) + a32*y2(-1) + a33*y3(-1)                   + b33*e3 ;
end;

histval;
    y1(0) = 0;
    y2(0) = 0;
    y3(0) = 0;
end;

shocks;
    var e1 = 1.0;
    var e2 = 1.0;
    var e3 = 1.0;
end;

steady;

check;

if isoctave 
    options_.bnlms.set_dynare_seed_to_default=false;
    set_dynare_seed(1);
end

TrueData = simul_backward_model([], 200);
// Set the periods where some of the endogenous variables will be constrained.
subsample = 3Y:100Y;

// Load the generated data
SimulatedData = copy(TrueData);

// Set the constrained paths for the endogenous variables (Output and PhysicalCapitalStock).
constrainedpaths = SimulatedData{'y1'}(subsample);

// Set the instruments (innovations used to control the paths for the endogenous variables).
exogenousvariables = dseries([NaN(100, 1) TrueData{'e2','e3'}.data(1:100,:)], '1Y', M_.exo_names);

// Invert the model by calling the model_inversion routine.
[endogenousvariables, exogenousvariables] = model_inversion(constrainedpaths, exogenousvariables, SimulatedData, M_, options_, oo_);

// Check that all the constraints are satisfied.
if max(abs(constrainedpaths(subsample).y1.data-endogenousvariables(subsample).y1.data))>1e-12
   error('Constraint on y1 path is not satisfied!')
end

if max(abs(exogenousvariables(subsample).e2.data-SimulatedData(subsample).e2.data))>1e-12
   error('Constraint on e2 path is not satisfied!')
end

if max(abs(exogenousvariables(subsample).e3.data-SimulatedData(subsample).e3.data))>1e-12
   error('Constraint on e3 path is not satisfied!')
end

// Check consistency of the results.
if max(abs(SimulatedData(subsample).y2.data-endogenousvariables(subsample).y2.data))>1e-12
  error('Model inversion is not consistent with respect to y2')
end

if max(abs(SimulatedData(subsample).y3.data-endogenousvariables(subsample).y3.data))>1e-12
   error('Model inversion is not consistent with respect to y3')
end

if max(abs(exogenousvariables(subsample).e1.data-SimulatedData(subsample).e1.data))>1e-12
   error('Model inversion is not consistent with true innovations (e1)')
end
