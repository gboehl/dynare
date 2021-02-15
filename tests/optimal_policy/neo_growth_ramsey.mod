/*
 * Copyright (C) 2021 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */
/*
 * This file computes a second-order approximation of the neo-classical growth model.
 * It assesses the conditional and unconditional welfares computed by the evaluate_planner_objective function
 * and compares them to a by-hand assessment stemming from the results the model neo_growth.mod incur.
 */

var k z c;

varexo e;

parameters beta gamma alpha delta rho s;

beta = 0.987;
gamma = 1;
delta = 0.012;
alpha = 0.4;
rho = 0.95;
s = 0.007;

model;
k=exp(z)*k(-1)^(alpha)-c+(1-delta)*k(-1);
z=rho*z(-1)+s*e;
end;

steady_state_model;
z = 0;
end;

shocks;
var e;
stderr 1;
end;

planner_objective ln(c);
ramsey_model(instruments=(k,c), planner_discount=beta);

initval;
   k = ((1/beta-(1-delta))/alpha)^(1/(alpha-1));
   c = k^alpha-delta*k;
end;

steady;
resid;
stoch_simul(order=2, irf=0);

planner_objective_value = evaluate_planner_objective(M_, options_, oo_);

if ~exist('neo_growth_results.mat','file');
   error('neo_growth must be run first');
end;

oo1 = load('neo_growth_results','oo_');
M1 = load('neo_growth_results','M_');
options1 = load('neo_growth_results','options_');
unc_W_hand = oo1.oo_.mean(strmatch('W',M1.M_.endo_names,'exact'));

initial_condition_states = repmat(oo1.oo_.dr.ys,1,M1.M_.maximum_lag);
shock_matrix = zeros(1,M1.M_.exo_nbr);
y_sim = simult_(M1.M_,options1.options_,initial_condition_states,oo1.oo_.dr,shock_matrix,options1.options_.order);
cond_W_hand=y_sim(strmatch('W',M1.M_.endo_names,'exact'),2);

if abs((unc_W_hand - planner_objective_value(1))/unc_W_hand) > 1e-6;
   error('Inaccurate unconditional welfare assessment');
end;
if abs((cond_W_hand - planner_objective_value(2))/cond_W_hand) > 1e-6;
   error('Inaccurate conditional welfare assessment');
end;
