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
 * This file computes a second-order approximation of the neo-classical growth model.
 * It is called by neo_growth_ramsey.mod to compare by-hand calculations of unconditional
 * and conditional welfares and the output of the evaluate_planner_objective function.
 */
@#include "neo_growth_common.inc"

shocks;
var e;
stderr 1;
end;

steady;
resid;

stoch_simul(order=2, irf=0);
