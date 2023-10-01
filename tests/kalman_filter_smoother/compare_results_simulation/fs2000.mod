/*
 * This file replicates the estimation of the cash in advance model described
 * Frank Schorfheide (2000): "Loss function-based evaluation of DSGE models",
 * Journal of Applied Econometrics, 15(6), 645-670.
 *
 * The data are in file "fsdat_simul.m", and have been artificially generated.
 * They are therefore different from the original dataset used by Schorfheide.
 *
 * The equations are taken from J. Nason and T. Cogley (1994): "Testing the
 * implications of long-run neutrality for monetary business cycle models",
 * Journal of Applied Econometrics, 9, S37-S70.
 * Note that there is an initial minus sign missing in equation (A1), p. S63.
 *
 * This implementation was written by Michel Juillard. Please note that the
 * following copyright notice only applies to this Dynare implementation of the
 * model.
 */

/*
 * Copyright © 2004-2019 Dynare Team
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

var m P c e W R k d n l gy_obs gp_obs y dA;
varexo e_a e_m;

parameters alp bet gam mst rho psi del;

alp = 0.33;
bet = 0.99;
gam = 0.003;
mst = 1.011;
rho = 0.7;
psi = 0.787;
del = 0.02;

model;
dA = exp(gam+e_a);
log(m) = (1-rho)*log(mst) + rho*log(m(-1))+e_m;
-P/(c(+1)*P(+1)*m)+bet*P(+1)*(alp*exp(-alp*(gam+log(e(+1))))*k^(alp-1)*n(+1)^(1-alp)+(1-del)*exp(-(gam+log(e(+1)))))/(c(+1)*P(+1)*m(+1))=0;
W = l/n;
-(psi/(1-psi))*(c*P/(1-n))+l/n = 0;
R = P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(-alp)/W;
1/(c*P)-bet*P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)/(m*l*c(+1)*P(+1)) = 0;
c+k = exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)+(1-del)*exp(-(gam+e_a))*k(-1);
P*c = m;
m-1+d = l;
e = exp(e_a);
y = k(-1)^alp*n^(1-alp)*exp(-alp*(gam+e_a));
exp(gy_obs) = dA*y/y(-1);
exp(gp_obs) = (P/P(-1))*m(-1)/dA;
end;

shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
end;

steady_state_model;
  dA = exp(gam);
  gst = 1/dA;
  m = mst;
  khst = ( (1-gst*bet*(1-del)) / (alp*gst^alp*bet) )^(1/(alp-1));
  xist = ( ((khst*gst)^alp - (1-gst*(1-del))*khst)/mst )^(-1);
  nust = psi*mst^2/( (1-alp)*(1-psi)*bet*gst^alp*khst^alp );
  n  = xist/(nust+xist);
  P  = xist + nust;
  k  = khst*n;

  l  = psi*mst*n/( (1-psi)*(1-n) );
  c  = mst/P;
  d  = l - mst + 1;
  y  = k^alp*n^(1-alp)*gst^alp;
  R  = mst/bet;
  W  = l/n;
  ist  = y-c;
  q  = 1 - d;

  e = 1;
  
  gp_obs = log(m/dA);
  gy_obs = log(dA);
end;

steady;

check;

estimated_params;
alp, beta_pdf, 0.356, 0.02;
bet, beta_pdf, 0.993, 0.002;
gam, normal_pdf, 0.0085, 0.003;
mst, normal_pdf, 1.0002, 0.007;
rho, beta_pdf, 0.129, 0.05;
psi, beta_pdf, 0.65, 0.05;
del, beta_pdf, 0.01, 0.005;
stderr e_a, inv_gamma_pdf, 0.035449, inf;
stderr e_m, inv_gamma_pdf, 0.008862, inf;
end;

varobs gp_obs gy_obs;

estimation(order=1,datafile=fsdat_simul_logged, silent_optimizer,consider_all_endogenous,nobs=192,mh_replic=2000, mh_nblocks=1,smoother, mh_jscale=0.8);

ex_=[];
for shock_iter=1:M_.exo_nbr
ex_=[ex_ oo_.SmoothedShocks.Mean.(M_.exo_names{shock_iter})];
end

ex_ = ex_(2:end,:);
% ex_ = zeros(size(ex_));
y0=[];
for endo_iter=1:M_.endo_nbr
y0 = [y0; oo_.SmoothedVariables.Mean.(M_.endo_names{endo_iter})(1)];
end;

%make sure decision rules were updated
[oo_.dr,info,M_.params] = resol(0,M_, options_, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);

dr = oo_.dr;
iorder=1;
y_=simult_(M_,options_,y0,dr,ex_,iorder);

fsdat_simul_logged;

%Needs bigger tolerance than ML, because transformation from parameters to steady states is not linear and steady state at mean parameters is not mean of steady states
if mean(abs(y_(strmatch('gy_obs',M_.endo_names,'exact'),:)'-(gy_obs(1:options_.nobs))))>1e-3 ||...
    mean(abs(y_(strmatch('gy_obs',M_.endo_names,'exact'),:)'-oo_.SmoothedVariables.Mean.gy_obs))>1e-3 ||...
    mean(abs(y_(strmatch('gp_obs',M_.endo_names,'exact'),:)'-(gp_obs(1:options_.nobs))))>1e-1 ||...
    mean(abs(y_(strmatch('gp_obs',M_.endo_names,'exact'),:)'-oo_.SmoothedVariables.Mean.gp_obs))>1e-2 
error('Smoother is wrong')
end


% figure
% plot((gy_obs))
% hold on
% plot(y_(strmatch('gy_obs',M_.endo_names,'exact'),:),'r--')
% 
% figure
% plot((gp_obs))
% hold on
% plot(y_(strmatch('gp_obs',M_.endo_names,'exact'),:),'r--')