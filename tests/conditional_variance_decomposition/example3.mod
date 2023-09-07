/*
 * Example 1 from F. Collard (2001): "Stochastic simulations with DYNARE:
 * A practical guide" (see "guide.pdf" in the documentation directory).
 * 
 * This file uses the steady_state_model-block to provide analytical steady state values.
 * To do so, the equations of the model have been transformed into a non-linear equation in 
 * labor h. Within the steady_state_model-block, a helper function is called that uses fsolve
 * to solve this non-linear equation. The use of the helper function is necessary to avoid 
 * interference of the MATLAB syntax with Dynare's preprocessor. A more complicated alternative 
 * that provides more flexibility in the type of commands executed and functions called is the use 
 * of an explicit steady state file. See the NK_baseline.mod in the Examples Folder.
 * 
 * This mod-file also shows how to use Dynare's capacities to generate TeX-files of the model equations. 
 * If you want to see the model equations belonging to this mod-file, run it using Dynare 
 * and then use a TeX-editor to compile the TeX-files generated.
 */

/*
 * Copyright © 2013 Dynare Team
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

@#define unit_root_var=0

var y, c, k, a, h, b
@#if unit_root_var==1
    , unit_root
@#endif
;
varexo e, u;

parameters beta $\beta$
     rho $\rho$
     alpha $\alpha$
     delta $\delta$
     theta $\theta$
     psi $\psi$
     tau $\tau$;

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
k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
    *(exp(b(+1))*alpha*y(+1)+(1-delta)*k));
y = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
k = exp(b)*(y-c)+(1-delta)*k(-1);
a = rho*a(-1)+tau*b(-1) + e;
b = tau*a(-1)+rho*b(-1) + u;
@#if unit_root_var==1
    unit_root=unit_root(-1)+e;
@#endif
end;

steady_state_model;
h=example3_steady_state_helper(alpha,beta,delta,psi,theta);
k=((1/beta-(1-delta))/alpha)^(1/(alpha-1))*h;
y = k^alpha*h^(1-alpha);
c=(1-alpha)*y/(theta*h^(1+psi));
a=0;
b=0;
end;

shocks;
var e; stderr 0.009;
var u; stderr 0.009;
var e, u = phi*0.009*0.009;
end;

stoch_simul(irf=0,conditional_variance_decomposition=[1,4,40],pruning,order=1);
oo1_=oo_;
stoch_simul(irf=0,conditional_variance_decomposition=[1,4,40],pruning,order=1) y k;
oo2_=oo_;
stoch_simul(irf=0,conditional_variance_decomposition=[1,4,40],pruning,order=2) y k ;
oo3_=oo_;
stoch_simul(irf=0,conditional_variance_decomposition=[1,4,40],pruning,order=2);
oo4_=oo_;

if max(max(abs(oo1_.variance_decomposition-oo4_.variance_decomposition)))>1e-8 || max(max(abs(oo2_.variance_decomposition-oo3_.variance_decomposition)))>1e-8
    error('Unconditional variance decomposition does not match.')
end

if max(max(max(abs(oo1_.conditional_variance_decomposition-oo4_.conditional_variance_decomposition))))>1e-8 || max(max(max(abs(oo2_.conditional_variance_decomposition-oo3_.conditional_variance_decomposition)))) >1e-8
    error('Conditional variance decomposition does not match.')
end

varobs y;
shocks;
var y; stderr 0.01;
end;

stoch_simul(irf=0,conditional_variance_decomposition=[1,4,40],pruning,order=1);
oo1_=oo_;
stoch_simul(irf=0,conditional_variance_decomposition=[1,4,40],pruning,order=1) y k;
oo2_=oo_;
stoch_simul(irf=0,conditional_variance_decomposition=[1,4,40],pruning,order=2) y k ;
oo3_=oo_;
stoch_simul(irf=0,conditional_variance_decomposition=[1,4,40],pruning,order=2);
oo4_=oo_;


if max(max(abs(oo1_.variance_decomposition-oo4_.variance_decomposition)))>1e-8 || max(max(abs(oo2_.variance_decomposition-oo3_.variance_decomposition)))>1e-8
    error('Unconditional variance decomposition does not match.')
end

if max(max(max(abs(oo1_.conditional_variance_decomposition-oo4_.conditional_variance_decomposition))))>1e-8 || max(max(max(abs(oo2_.conditional_variance_decomposition-oo3_.conditional_variance_decomposition)))) >1e-8
    error('Conditional variance decomposition does not match.')
end

if max(max(abs(oo1_.variance_decomposition_ME-oo4_.variance_decomposition_ME)))>1e-2 || max(max(abs(oo2_.variance_decomposition_ME-oo3_.variance_decomposition_ME)))>1e-2
    error('Unconditional variance decomposition with ME does not match.')
end

if max(max(max(abs(oo1_.conditional_variance_decomposition_ME-oo4_.conditional_variance_decomposition_ME))))>1e-8 || max(max(max(abs(oo2_.conditional_variance_decomposition_ME-oo3_.conditional_variance_decomposition_ME))))>1e-8
    error('Conditional variance decomposition with ME does not match.')
end
