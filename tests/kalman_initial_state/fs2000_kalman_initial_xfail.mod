@#include "fs2000_common.inc"

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

options_.solve_tolf = 1e-12;

filter_initial_state;
k(0)= (((1-1/exp(gam)*bet*(1-del)) / (alp*1/exp(gam)^alp*bet) )^(1/(alp-1)))*(( (((((1-1/exp(gam)*bet*(1-del)) / (alp*1/exp(gam)^alp*bet) )^(1/(alp-1)))*1/exp(gam))^alp - (1-1/exp(gam)*(1-del))*(((1-1/exp(gam)*bet*(1-del)) / (alp*1/exp(gam)^alp*bet) )^(1/(alp-1))))/mst )^(-1))/(psi*mst^2/( (1-alp)*(1-psi)*bet*1/exp(gam)^alp*(((1-1/exp(gam)*bet*(1-del)) / (alp*1/exp(gam)^alp*bet) )^(1/(alp-1)))^alp )+(( (((((1-1/exp(gam)*bet*(1-del)) / (alp*1/exp(gam)^alp*bet) )^(1/(alp-1)))*1/exp(gam))^alp - (1-1/exp(gam)*(1-del))*(((1-1/exp(gam)*bet*(1-del)) / (alp*1/exp(gam)^alp*bet) )^(1/(alp-1))))/mst)^(-1)));
c(0)=1;
end;

estimation(order=1,datafile='../fs2000/fsdat_simul',nobs=192,loglinear,mh_replic=2000,mh_nblocks=1,mh_jscale=0.8,moments_varendo,consider_only_observed);
