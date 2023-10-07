/*
 * This file is based on the cash in advance model described
 * Frank Schorfheide (2000): "Loss function-based evaluation of DSGE models",
 * Journal of Applied Econometrics, 15(6), 645-670.
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
 * Copyright © 2004-2013 Dynare Team
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

parameters alp bet gam mst rho psi del theta;

alp = 0.33;
bet = 0.99;
gam = 0.003;
mst = 1.011;
rho = 0.7;
psi = 0.787;
del = 0.02;
theta=0;

model;
dA = exp(gam+e_a);
log(m) = (1-rho)*log(mst) + rho*log(m(-1))+e_m;
-P/(c(+1)*P(+1)*m)+bet*P(+1)*(alp*exp(-alp*(gam+log(e(+1))))*k^(alp-1)*n(+1)^(1-alp)+(1-del)*exp(-(gam+log(e(+1)))))/(c(+2)*P(+2)*m(+1))=0;
W = l/n;
-(psi/(1-psi))*(c*P/(1-n))+l/n = 0;
R = P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(-alp)/W;
1/(c*P)-bet*P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)/(m*l*c(+1)*P(+1)) = 0;
c+k = exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)+(1-del)*exp(-(gam+e_a))*k(-1);
P*c = m;
m-1+d = l;
e = exp(e_a);
y = k(-1)^alp*n^(1-alp)*exp(-alp*(gam+e_a));
gy_obs = dA*y/y(-1);
gp_obs = (P/P(-1))*m(-1)/dA;
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
  
  gp_obs = m/dA;
  gy_obs = dA;
end;

varobs gp_obs gy_obs;

shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
corr gy_obs,gp_obs = 0.25;
end;

stoch_simul(order=1,irf=0,nomoments,nofunctions,periods=1000);

verbatim;
temp=oo_.endo_simul([strmatch('gy_obs',M_.endo_names,'exact'),strmatch('gp_obs',M_.endo_names,'exact')],:);
exo_simul = randn(2,options_.periods)'*chol([0.005 0.25*sqrt(0.005*0.01); 0.25*sqrt(0.005*0.01) 0.01]);
temp=temp'+exo_simul;
gy_obs=temp(:,1);
gp_obs=temp(:,2);
save('fsdat_simul_corr.mat','gy_obs','gp_obs')
end;
estimated_params;
alp, 0.33, 0.2, 0.4;
gam,  0.0085, 0 , 0.1;
del, 0.01, 0, 0.1;
stderr e_a, 0.014, 0, Inf;
stderr e_m, 0.005, 0, Inf;
corr e_m, e_a, 0, 0, 1;
stderr gp_obs, 1, 0, Inf;
stderr gy_obs, 1, 0, Inf;
corr gp_obs, gy_obs,0, 0 , 1;
end;
warning('off','MATLAB:nearlySingularMatrix')
estimation(order=1,datafile=fsdat_simul_corr,silent_optimizer,prior_trunc=0,mode_check,smoother,filter_decomposition,forecast = 8,filtered_vars,filter_step_ahead=[1,3],irf=20,tex) m P c e W R k d y gy_obs;



estimated_params(overwrite);
//alp, beta_pdf, 0.356, 0.02;
gam, normal_pdf, 0.0085, 0.003;
//del, beta_pdf, 0.01, 0.005;
stderr e_a, inv_gamma_pdf, 0.05449, inf;
stderr e_m, inv_gamma_pdf, 0.008862, inf;
corr e_m, e_a, , -1, 1, normal_pdf, 0, 0.2;
stderr gp_obs, 0.1, inv_gamma_pdf, 0.001, inf;
//stderr gy_obs, inv_gamma_pdf, 0.001, inf;
//corr gp_obs, gy_obs,normal_pdf, 0, 0.2;
end;

estimation(mode_compute=5,silent_optimizer,order=1,datafile=fsdat_simul_corr,mode_check,smoother,filter_decomposition,mh_replic=2000, mh_nblocks=1, mh_jscale=0.8,forecast = 8,bayesian_irf,filtered_vars,filter_step_ahead=[1,3],irf=20) m P c e W R k d y;
shock_decomposition y W R;
//identification(advanced=1,max_dim_cova_group=3,prior_mc=250);
