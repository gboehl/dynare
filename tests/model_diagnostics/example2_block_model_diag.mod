//test singular Jacobian check with block decomposition

/*
 * Example 2 from F. Collard (2001): "Stochastic simulations with DYNARE:
 * A practical guide" (see "guide.pdf" in the documentation directory).
 */

/*
 * Copyright (C) 2001-2010 Dynare Team
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


var y, c, k, a, h, b;
varexo e, u;

parameters beta, rho, alpha, delta, theta, psi, tau;

alpha = 0.36;
rho   = 1;
tau   = 0;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;

model(block);
exp(c)*theta*exp(h)^(1+psi)=(1-alpha)*exp(y);
exp(k) = beta*(((exp(b)*exp(c))/(exp(b(+1))*exp(c(+1))))
         *(exp(b(+1))*alpha*exp(y(+1))+(1-delta)*exp(k)));
exp(y) = exp(a)*(exp(k(-1))^alpha)*(exp(h)^(1-alpha));
exp(k) = exp(b)*(exp(y)-exp(c))+(1-delta)*exp(k(-1));
a = rho*a(-1)+tau*b(-1) + e;
b = tau*a(-1)+rho*b(-1) + u;
end;

initval;
y = 0.1;
c = -0.2;
h = -1.2;
k =  2.4;
a = 0;
b = 0;
e = 0;
u = 0;
end;

steady;
model_diagnostics;