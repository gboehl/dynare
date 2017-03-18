//conducts posterior sampling with morris=1 and morris=2 and rmse
var y l w pi r A dy pic;
varexo eps_A eps_P;

parameters PSI ETA BETA ALPHA KAPPA PI RHO;

PSI = 1;
ETA = 2;
BETA = 0.95;
KAPPA =  0.5;
ALPHA =  0.99;
PI = 1.005;
RHO = 0.9;

model;
y = exp(A)*l; //production function
w/y = PSI * l^ETA;//labor supply
BETA*1/y(+1)*r/pi(+1)=1/y;//saving
log(pi) =  ALPHA*log(pi(-1))+(1-ALPHA) * log(pi(+1)) + KAPPA*log(w/exp(A))+eps_P;//inflation
log(r*BETA/PI) = 1.5*log(pi/PI);
A = RHO*A(-1)+eps_A;
//observable variables
dy = log(y)-log(y(-1));
pic= log(pi)-log(PI);
end;


steady_state_model;
A = eps_A;
w = exp(A);
pi = PI;
l = (1/PSI)^(1/(1+ETA));
y = exp(A)*l;
r = PI/BETA;
end;

shocks;
var eps_A;
stderr 0.01;
var eps_P;
stderr 0.01;

end;

options_.TeX=1;

steady;
check;

varobs dy pic;
estimated_params; 
ETA,3.7,0.00000001,10,GAMMA_PDF,5,1;
KAPPA,0.3,0.00000001,0.99999999999,beta_PDF,0.5,0.1;
stderr eps_A,0.02,0.000000000001,100,INV_GAMMA2_PDF,0.2,inf;
stderr eps_P,0.03,0.000000000001,100,INV_GAMMA2_PDF,0.2,inf;
end;
estimation(order=1,prior_trunc=0,plot_priors =0, datafile=nk_est_data,conf_sig =.95,smoother,moments_varendo,filtered_vars,mode_check,mode_compute=4,mh_replic=5000,mh_jscale=1.5,mh_nblocks=1,bayesian_irf) y pi l dy pic;
dynare_sensitivity (datafile=nk_est_data,rmse=0, nsam = 2000, lik_only = 0, morris = 2,var_rmse=(dy pic)) ;
dynare_sensitivity (datafile=nk_est_data,rmse=0, nsam = 2000, lik_only = 0, morris = 1,var_rmse=(dy pic)) ;


