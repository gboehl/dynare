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

var y1 y2 y3 ;

varexo e1 e2 e3 ;

parameters a11 a12 a13 a21 a22 a23 a31 a32 a33 ;

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

model(linear);
    diff(y1) = a11*diff(y1(-1)) + a12*diff(y2(-1)) + a13*diff(y3(-1)) + e1 ;
    diff(y2) = a21*diff(y1(-1)) + a22*diff(y2(-1)) + a23*diff(y3(-1)) + e2 ;
    diff(y3) = a31*diff(y1(-1)) + a32*diff(y2(-1)) + a33*diff(y3(-1)) + e3 ;
end;

shocks;
    var e1 = 1.0;
    var e2 = 1.0;
    var e3 = 1.0;
end;

steady;

check;

@#if USE_HISTVAL
    histval;
        y1(0) = 0;
    	y2(0) = 0;
    	y3(0) = 0;
    	y1(-1) = 0;
    	y2(-1) = 0;
    	y3(-1) = 0;
    end;
    TrueData = simul_backward_model([], 200);
@#else
    initialconditions = dseries(zeros(2,6), 1Y, vertcat(M_.endo_names(1:M_.orig_endo_nbr), M_.exo_names));
    TrueData = simul_backward_model(initialconditions, 200);
@#endif

// Set the periods where some of the endogenous variables will be constrained.
subsample = 3Y:100Y;

// Load the generated data
SimulatedData = copy(TrueData);

// Set the constrained paths for the endogenous variables (Output and PhysicalCapitalStock).
constrainedpaths = SimulatedData{'y1'}(subsample);

// Set the instruments (innovations used to control the paths for the endogenous variables).
exogenousvariables = dseries([NaN(98, 1) TrueData{'e2','e3'}.data(3:100,:)], 3Y, M_.exo_names);

// Invert the model by calling the model_inversion routine.
[endogenousvariables, exogenousvariables] = model_inversion(constrainedpaths, exogenousvariables, SimulatedData, M_, options_, oo_);

// Check that all the constraints are satisfied.
if max(abs(constrainedpaths(subsample).y1.data-endogenousvariables(subsample).y1.data))>1e-12
   error('Constraint on y1 path is not satisfied!')
end

if max(abs(exogenousvariables(subsample).e2.data-SimulatedData(subsample).e2.data))>1e-12
   error('Constraint on e1 path is not satisfied!')
end

if max(abs(exogenousvariables(subsample).e3.data-SimulatedData(subsample).e3.data))>1e-12
   error('Constraint on e2 path is not satisfied!')
end

// Check consistency of the results.
if max(abs(SimulatedData(subsample).y2.data-endogenousvariables(subsample).y2.data))>1e-12
  error('Model inversion is not consistent with respect to y2')
end

if max(abs(SimulatedData(subsample).y3.data-endogenousvariables(subsample).y3.data))>1e-12
   error('Model inversion is not consistent with respect to y3')
end

if max(abs(exogenousvariables(subsample).e1.data-SimulatedData(subsample).e1.data))>1e-12
   error('Model inversion is not consistent with true innovations (e3)')
end
