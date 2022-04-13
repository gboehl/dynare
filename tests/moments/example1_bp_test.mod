/*
 * Example 1 from F. Collard (2001): "Stochastic simulations with DYNARE:
 * A practical guide" (see "guide.pdf" in the documentation directory).
 */

/*
 * Copyright © 2001-2020 Dynare Team
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
rho   = 0.5;
tau   = 0.025;
beta  = 0.98;
delta = 0.025;
psi   = 0;
theta = 2.95;

phi   = 0.1;

model;
c*theta*h^(1+psi)=(1-alpha)*y;
k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
    *(exp(b(+1))*alpha*y(+1)+(1-delta)*k));
y = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
k = exp(b)*(y-c)+(1-delta)*k(-1);
a = rho*a(-1)+tau*b(-1) + e;
b = tau*a(-1)+rho*b(-1) + u;
end;

initval;
y = 1.08068253095672;
c = 0.80359242014163;
h = 0.29175631001732;
k = 11.08360443260358;
a = 0;
b = 0;
e = 0;
u = 0;
end;

shocks;
var e; stderr 0.009;
var u; stderr 0.009;
var e, u = phi*0.009*0.009;
end;

steady(solve_algo=4,maxit=1000);

stoch_simul(order=1,nofunctions,irf=0,bandpass_filter=[6 32],filtered_theoretical_moments_grid=8192);
oo_filtered_all_shocks_theoretical=oo_;

stoch_simul(order=1,nofunctions,periods=1000000);
oo_filtered_all_shocks_simulated=oo_;


shocks;
var e; stderr 0;
var u; stderr 0.009;
var e, u = phi*0.009*0;
end;

stoch_simul(order=1,nofunctions,irf=0,periods=0);

oo_filtered_one_shock_theoretical=oo_;

stoch_simul(order=1,nofunctions,periods=1000000);
oo_filtered_one_shock_simulated=oo_;


if max(abs((1-diag(oo_filtered_one_shock_simulated.var)./(diag(oo_filtered_all_shocks_simulated.var)))*100-oo_filtered_all_shocks_theoretical.variance_decomposition(:,1)))>2
    error('Variance Decomposition wrong')
end

if max(max(abs(oo_filtered_all_shocks_simulated.var-oo_filtered_all_shocks_theoretical.var)))>5e-4;
    error('Covariance wrong')
end

if max(max(abs(oo_filtered_one_shock_simulated.var-oo_filtered_one_shock_theoretical.var)))>5e-4;
    error('Covariance wrong')
end

for ii=1:options_.ar
    autocorr_model_all_shocks_simulated(:,ii)=diag(oo_filtered_all_shocks_simulated.autocorr{ii});
    autocorr_model_all_shocks_theoretical(:,ii)=diag(oo_filtered_all_shocks_theoretical.autocorr{ii});
    autocorr_model_one_shock_simulated(:,ii)=diag(oo_filtered_one_shock_simulated.autocorr{ii});
    autocorr_model_one_shock_theoretical(:,ii)=diag(oo_filtered_one_shock_theoretical.autocorr{ii});
end

if max(max(abs(autocorr_model_all_shocks_simulated-autocorr_model_all_shocks_theoretical)))>0.05;
    error('Correlation wrong')
end

if max(max(abs(autocorr_model_one_shock_simulated-autocorr_model_one_shock_theoretical)))>0.05;
    error('Correlation wrong')
end
