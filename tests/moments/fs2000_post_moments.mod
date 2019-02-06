/*
 * This file replicates the estimation of the cash in advance model (termed M1 
 * in the paper) described in Frank Schorfheide (2000): "Loss function-based 
 * evaluation of DSGE models", Journal of Applied Econometrics, 15(6), 645-670.
 *
 * The data are in file "fsdat_simul.m", and have been artificially generated.
 * They are therefore different from the original dataset used by Schorfheide.
 *
 * The prior distribution follows the one originally specified in Schorfheide's
 * paper, except for parameter rho. In the paper, the elicited beta prior for rho
 * implies an asymptote and corresponding prior mode at 0. It is generally
 * recommended to avoid this extreme type of prior. Some optimizers, for instance
 * mode_compute=12 (Mathworks' particleswarm algorithm) may find a posterior mode
 * with rho equal to zero. We lowered the value of the prior standard deviation
 * (changing .223 to .100) to remove the asymptote.
 *
 * The equations are taken from J. Nason and T. Cogley (1994): "Testing the
 * implications of long-run neutrality for monetary business cycle models",
 * Journal of Applied Econometrics, 9, S37-S70.
 * Note that there is an initial minus sign missing in equation (A1), p. S63.
 *
 * This implementation was originally written by Michel Juillard. Please note that the
 * following copyright notice only applies to this Dynare implementation of the
 * model.
 */

/*
 * Copyright (C) 2004-2017 Dynare Team
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
  
  gp_obs = m/dA;
  gy_obs = dA;
end;

steady;

check;

estimated_params;
alp, beta_pdf, 0.356, 0.02;
bet, beta_pdf, 0.993, 0.002;
gam, normal_pdf, 0.0085, 0.003;
mst, normal_pdf, 1.0002, 0.007;
rho, beta_pdf, 0.129, 0.100;
psi, beta_pdf, 0.65, 0.05;
del, beta_pdf, 0.01, 0.005;
stderr e_a, inv_gamma_pdf, 0.035449, inf;
stderr e_m, inv_gamma_pdf, 0.008862, inf;
end;

varobs gp_obs gy_obs;

estimation(order=1,mode_compute=5, datafile='../fs2000/fsdat_simul.m', nobs=192, loglinear, mh_replic=20, mh_nblocks=1, mh_jscale=0.8,moments_varendo,
conditional_variance_decomposition=[2,2000],consider_all_endogenous,sub_draws=2);

stoch_simul(order=1,conditional_variance_decomposition=[2,2000],noprint,nograph);
par=load([M_.fname filesep 'metropolis' filesep M_.fname '_posterior_draws1']);

for par_iter=1:size(par.pdraws,1)
   M_=set_parameters_locally(M_,par.pdraws{par_iter,1});
   info=stoch_simul(var_list_);
   correlation(:,:,par_iter)=cell2mat(oo_.autocorr);
   covariance(:,:,par_iter)=oo_.var;
   conditional_variance_decomposition(:,:,:,par_iter)=oo_.conditional_variance_decomposition;
   variance_decomposition(:,:,par_iter)=oo_.variance_decomposition;
end

correlation=mean(correlation,3);
nvars=M_.orig_endo_nbr;
for var_iter_1=1:nvars
    for var_iter_2=1:nvars
        if max(abs(correlation(var_iter_1,var_iter_2:nvars:end)'-oo_.PosteriorTheoreticalMoments.dsge.correlation.Mean.(deblank(M_.endo_names(var_iter_1,:))).(deblank(M_.endo_names(var_iter_2,:)))))>1e-8
            error('Correlations do not match')
        end
    end
end

covariance=mean(covariance,3);
nvars=size(M_.endo_names(1:M_.orig_endo_nbr,:),1);
for var_iter_1=1:nvars
    for var_iter_2=var_iter_1:nvars
        if max(abs(covariance(var_iter_1,var_iter_2)-oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.(deblank(M_.endo_names(var_iter_1,:))).(deblank(M_.endo_names(var_iter_2,:)))))>1e-8
            error('Covariances do not match')
        end
    end
end

variance_decomposition=mean(variance_decomposition,3);
nvars=size(M_.endo_names(1:M_.orig_endo_nbr,:),1);
for var_iter_1=1:nvars
    for shock_iter=1:M_.exo_nbr
        if max(abs(variance_decomposition(var_iter_1,shock_iter)/100-oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.(deblank(M_.endo_names(var_iter_1,:))).(deblank(M_.exo_names(shock_iter,:)))))>1e-8
            error('Variance decomposition does not match')
        end
    end
end

conditional_variance_decomposition=mean(conditional_variance_decomposition,4);
nvars=size(M_.endo_names(1:M_.orig_endo_nbr,:),1);
horizon_size=size(conditional_variance_decomposition,3);
for var_iter_1=1:nvars
    for shock_iter=1:M_.exo_nbr
        for horizon_iter=1:horizon_size
            if max(abs(conditional_variance_decomposition(var_iter_1,horizon_iter,shock_iter)-oo_.PosteriorTheoreticalMoments.dsge.ConditionalVarianceDecomposition.Mean.(deblank(M_.endo_names(var_iter_1,:))).(deblank(M_.exo_names(shock_iter,:)))(horizon_iter)))>1e-8
                error('Conditional Variance decomposition does not match')
            end
        end
    end
end

/*
 * The following lines were used to generate the data file. If you want to
 * generate another random data file, comment the "estimation" line and uncomment
 * the following lines.
 */

//stoch_simul(periods=200, order=1);
//datatomfile('fsdat_simul', char('gy_obs', 'gp_obs'));
