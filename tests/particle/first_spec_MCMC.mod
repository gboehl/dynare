//File testing error message if initial state vector is not positive definite

@#include "first_spec_common.inc"

varobs b ca;

shocks;
var eeps = 0.04^2;
var nnu = 0.03^2;
var b = 0.01^2;
var ca = 0.01^2;
end;

stoch_simul(order=3,periods=200, irf=0);

save('my_data_MCMC.mat','ca','b');

options_.pruning=1;
options_.particle.pruning=1;
options_.particle.number_of_particles=500;

estimation(datafile='my_data_MCMC.mat',order=2,mh_replic=100,filter_algorithm=sis,nonlinear_filter_initialization=2
    ,mode_compute=0 %don't compute mode
    ,mcmc_jumping_covariance=identity_matrix %use identity matrix
    ,cova_compute=0 %tell program that no covariance matrix was computed
    ,bayesian_irf,moments_varendo,consider_all_endogenous);