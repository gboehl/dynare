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
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */
/*
 * This file computes a kth-order approximation of the neo-classical growth model.
 * It assesses the conditional welfare and the derivatives of the felicity and welfare functions computed by the k_order_welfare function
 * and compares them to a by-hand assessment stemming from the results the model neo_growth_k_order.mod incur.
 */

@#include "neo_growth_ramsey_common.inc"

shocks;
var e;
stderr 1;
end;

stoch_simul(order=6, irf=0);

evaluate_planner_objective;

[W_dynpp] = k_order_welfare(oo_.dr, M_, options_);

if ~exist(['neo_growth_k_order' filesep 'Output' filesep 'neo_growth_k_order_results.mat'],'file');
   error('neo_growth_k_order must be run first');
end;

oo = load(['neo_growth_k_order' filesep 'Output' filesep 'neo_growth_k_order_results'],'oo_');
M = load(['neo_growth_k_order' filesep 'Output' filesep 'neo_growth_k_order_results'],'M_');
options = load(['neo_growth_k_order' filesep 'Output' filesep 'neo_growth_k_order_results'],'options_');

ind_W = strmatch('W', M.M_.endo_names,'exact');

err = -1e6;
for i = 1:options_.order
   field_W = strcat('W_', num2str(i));
   g_i = eval(strcat('oo.oo_.dr.g_', num2str(i)));
   tmp_err = max(W_dynpp.(field_W) - g_i(ind_W, :));
   err = max(err, tmp_err);
end

if err > 1e-10;
   error('Inaccurate assessment of the derivatives of the welfare function');
end;
