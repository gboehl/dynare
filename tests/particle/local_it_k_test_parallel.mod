/*
  Tests that the parallelized version of local_state_space_iteration_k and local_state_space_iteration_fortran
  (for k=3) return the same results.

  This file must be run after first_spec.mod (both are based on the same model).
*/

@#include "first_spec_common.inc"

varobs q ca;

shocks;
var eeps = 0.04^2;
var nnu = 0.03^2;
var q = 0.01^2;
var ca = 0.01^2;
end;

// Initialize various structures
estimation(datafile='my_data.mat',order=2,mode_compute=0,mh_replic=0,filter_algorithm=sis,nonlinear_filter_initialization=2
    ,cova_compute=0 %tell program that no covariance matrix was computed
);

stoch_simul(order=2, periods=200, irf=0, k_order_solver);

// Really perform the test

nparticles = options_.particle.number_of_particles;
nsims = 1e6/nparticles;

/* We generate particles using realistic distributions (though this is not
   strictly needed) */
state_idx = oo_.dr.order_var((M_.nstatic+1):(M_.nstatic+M_.npred+M_.nboth));
yhat = chol(oo_.var(state_idx,state_idx))*randn(M_.npred+M_.nboth, nparticles);
epsilon = chol(M_.Sigma_e)*randn(M_.exo_nbr, nparticles);

stoch_simul(order=3, periods=200, irf=0, k_order_solver);

options_.threads.local_state_space_iteration_k = 2;
options_.threads.local_state_space_iteration_fortran = 2;

tStart1 = tic; for i=1:nsims, ynext1 = local_state_space_iteration_k(yhat, epsilon, oo_.dr, M_, options_); end, tElapsed1 = toc(tStart1);

tStart2 = tic; [udr] = folded_to_unfolded_dr(oo_.dr, M_, options_); for i=1:nsims, ynext2 = local_state_space_iteration_fortran(yhat, epsilon, oo_.dr, M_, options_, udr); end, tElapsed2 = toc(tStart2);

if max(max(abs(ynext1-ynext2))) > 1e-14
    error('Inconsistency between local_state_space_iteration_k and local_state_space_iteration_fortran')
end

if tElapsed2<tElapsed1
    skipline()
    dprintf('local_state_space_iteration_fortran is %5.2f times faster than local_state_space_iteration_k', tElapsed1/tElapsed2)
    skipline()
else
    skipline()
    dprintf('local_state_space_iteration_fortran is %5.2f times slower than local_state_space_iteration_k', tElapsed2/tElapsed1)
    skipline()
end
