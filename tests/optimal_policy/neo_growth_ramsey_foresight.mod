/*
 * Copyright © 2021 Dynare Team
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
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */
 /*
 * This file simulates a perfect-foresight version of the neo-classical growth model.
 * It assesses the conditional and unconditional welfares computed by the evaluate_planner_objective function
 * and compares them to a by-hand assessment stemming from the results of the model neo_growth_foresight.mod
 */

@#include "neo_growth_ramsey_common.inc"

shocks;
var e;
periods 1;
values 1;
end;

perfect_foresight_setup(periods=200);
perfect_foresight_solver;

planner_objective_value = evaluate_planner_objective(M_, options_, oo_);

if ~exist('neo_growth_foresight_results.mat','file');
   error('neo_growth_foresight must be run first');
end;

oo1 = load(['neo_growth_foresight' filesep 'Output' filesep 'neo_growth_foresight_results'],'oo_');
M1 = load(['neo_growth_foresight' filesep 'Output' filesep 'neo_growth_foresight_results'],'M_');
options1 = load(['neo_growth_foresight' filesep 'Output' filesep 'neo_growth_foresight_results'],'options_');
cond_W_hand = oo1.oo_.endo_simul(strmatch('W',M1.M_.endo_names,'exact'),2);
unc_W_hand = oo1.oo_.endo_simul(strmatch('W',M1.M_.endo_names,'exact'),end);

if abs((unc_W_hand - planner_objective_value.unconditional)/unc_W_hand) > 1e-6;
   error('Inaccurate unconditional welfare assessment');
end;
if abs((cond_W_hand - planner_objective_value.conditional)/cond_W_hand) > 1e-6;
   error('Inaccurate conditional welfare assessment');
end;
