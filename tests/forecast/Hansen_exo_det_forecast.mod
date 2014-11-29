/*
* The Hansen model following McCandless, George T. (2008): The ABCs of RBCs, Hardvard University Press, Chapter 6
*
* This mod-file tests the correctness of forecasting with exogenous deterministic variables.
* A forecast starting at the steady state in t=0, where the only shock is a perfectly 
* anticipated shock in t=8, is equal to the IRF to a 7 period anticipated news shock.
* Note the timing difference due to the fact that in forecasting. the agent starts at the 
* steady state at time 0 and has the first endogenous reaction period at t=1 so that the shock at 
* t=8 is only 7 period anticipated

* This implementation was written by Johannes Pfeifer. Please note that the
* following copyright notice only applies to this Dynare implementation of the
* model.
*/

/*
 * Copyright (C) 2014 Johannes Pfeifer
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



var c r y h k a;
varexo eps_a_surp eps_a_antic;
varexo_det e_det;

parameters beta delta theta rho_a A_disutil;

beta = 0.99;
delta = 0.025;
theta = 0.36;
rho_a = 0.35;
A_disutil = 2;

model;
1/c = beta*((1/c(+1))*(r(+1) +(1-delta)));
(1-h)*(1-theta)*(y/h) = A_disutil*c;
c = y +(1-delta)*k(-1) - k;
y = a*k(-1)^(theta)*h^(1-theta);
r = theta*(y/k(-1));
log(a)=rho_a*log(a(-1))+eps_a_surp+eps_a_antic(-7)+e_det;
end;

steady_state_model;
a = 1;
h = (1+(A_disutil/(1-theta))*(1 - (beta*delta*theta)/(1-beta*(1-delta))))^(-1);
k = h*(theta*a/(1/beta -(1-delta)))^(1/(1-theta));
y = a*k^(theta)*h^(1-theta);
c = y-delta*k;
r =  1/beta - (1-delta);
end;

steady;

shocks;
var eps_a_surp; stderr 0;
var eps_a_antic; stderr 0.01;
var e_det;
periods 8;
values 0.01;
end;

stoch_simul(irf=40,nomoments, order=1);

% make_ex_;

forecast(periods=40);
figure
for ii=1:M_.orig_endo_nbr
    subplot(3,3,ii)    
    var_name=deblank(M_.endo_names(ii,:));
    var_index=strmatch(var_name,deblank(M_.endo_names),'exact');
    plot(1:40,oo_.dr.ys(var_index)+oo_.irfs.([var_name,'_eps_a_antic']),'b',1:40,... %
            oo_.forecast.Mean.(var_name),'r--')
    title(var_name)
    difference(ii,:)=oo_.dr.ys(var_index)+oo_.irfs.([var_name,'_eps_a_antic'])-oo_.forecast.Mean.(var_name)';
end

if max(max(abs(difference)))>1e-10
   error('Forecasts with exogenous deterministic variable is wrong') 
end