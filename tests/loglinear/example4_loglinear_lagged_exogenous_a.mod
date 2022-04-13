/*
 * Example 1 from F. Collard (2001): "Stochastic simulations with DYNARE:
 * A practical guide" (see "guide.pdf" in the documentation directory).
 */

/*
 * Copyright © 2001-2016 Dynare Team
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


var y, c, k, a, h, b, e1, u1;
varexo e, u;

parameters beta, rho, alpha, delta, theta, psi, tau;

alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;

phi   = 0.1;

model;
c*theta*h^(1+psi)=(1-alpha)*y;
k = beta*(((b*c)/(b(+1)*c(+1)))
    *(b(+1)*alpha*y(+1)+(1-delta)*k));
y = a*(k(-1)^alpha)*(h^(1-alpha));
k = b*(y-c)+(1-delta)*k(-1);
log(a) = rho*log(a(-1))+tau*log(b(-1)) + log(e1(-1));
log(b) = tau*log(a(-1))+rho*log(b(-1)) + log(u1(-1));
e1 = exp(e);
u1 = exp(u);
end;

initval;
y = 1.08068253095672;
c = 0.80359242014163;
h = 0.29175631001732;
k = 11.08360443260358;
a = 1;
b = 1;
e1 = 1;
u1 = 1;
end;

shocks;
var e; stderr 0.009;
var u; stderr 0.009;
var e, u = phi*0.009*0.009;
end;

stoch_simul(loglinear,order=1);

D = load(['example4_loglinear_lagged_exogenous' filesep 'Output' filesep 'example4_loglinear_lagged_exogenous_results']);

test1 = D.oo_.dr.ghx - oo_.dr.ghx;
if norm(test1) > 1e-16;
   error('error in computing ghx');
end;
test2 = D.oo_.dr.ghu - oo_.dr.ghu;
if norm(test2) > 1e-16;
   error('error in computing ghu');
end;

for i = fieldnames(D.oo_.irfs)';
    test3 = D.oo_.irfs.(i{1}) - oo_.irfs.(i{1});
    if norm(test2) > 1e-16;
        error(['error in computing irf ' i]);
    end;
end;        
    
