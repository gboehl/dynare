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
 * It is called by neo_growth_ramsey.mod to compare by-hand calculations of unconditional
 * and conditional welfares and the output of the evaluate_planner_objective function.
 */
var U k z c W;

varexo e;

parameters beta gamma alpha delta rho s;

beta = 0.987;
gamma = 1;
delta = 0.012;
alpha = 0.4;
rho = 0.95;
s = 0.007;

model;
c^(-gamma)=beta*c(+1)^(-gamma)*(alpha*exp(z(+1))*k^(alpha-1)+1-delta);
W = U + beta*W(+1);
k=exp(z)*k(-1)^(alpha)-c+(1-delta)*k(-1);
z=rho*z(-1)+s*e;
U=ln(c);
end;

steady_state_model;
k = ((1/beta-(1-delta))/alpha)^(1/(alpha-1));
c = k^alpha-delta*k;
z = 0;
U = ln(c);
W = U/(1-beta);
end;

shocks;
var e;
stderr 1;
end;

steady;
resid;
stoch_simul(order=2, irf=0);
