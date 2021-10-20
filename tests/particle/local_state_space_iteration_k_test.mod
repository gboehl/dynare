/*
  Tests that local_state_space_iteration_2 and local_state_space_iteration_k
  (for k=2) return the same results.

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

nparticles = 100;

/* We generate particles using realistic distributions (though this is not
   strictly needed) */
state_idx = oo_.dr.order_var((M_.nstatic+1):(M_.nstatic+M_.npred+M_.nboth));
yhat = chol(oo_.var(state_idx,state_idx))*randn(M_.npred+M_.nboth, nparticles);
epsilon = chol(M_.Sigma_e)*randn(M_.exo_nbr, nparticles);

dr = oo_.dr;

// “rf” stands for “Reduced Form”
rf_ghx = dr.ghx(dr.restrict_var_list, :);
rf_ghu = dr.ghu(dr.restrict_var_list, :);
rf_constant = dr.ys(dr.order_var)+0.5*dr.ghs2;
rf_constant = rf_constant(dr.restrict_var_list, :);
rf_ghxx = dr.ghxx(dr.restrict_var_list, :);
rf_ghuu = dr.ghuu(dr.restrict_var_list, :);
rf_ghxu = dr.ghxu(dr.restrict_var_list, :);

tStart1 = tic; for i=1:10000, ynext1 = local_state_space_iteration_2(yhat, epsilon, rf_ghx, rf_ghu, rf_constant, rf_ghxx, rf_ghuu, rf_ghxu, 1); end, tElapsed1 = toc(tStart1);

tStart2 = tic; for i=1:10000, ynext2 = local_state_space_iteration_k(yhat, epsilon, dr, M_, options_); end, tElapsed2 = toc(tStart2);


tStart3 = tic; [udr] = folded_to_unfolded_dr(dr, M_, options_); for i=1:10000, ynext3 = local_state_space_iteration_fortran(yhat, epsilon, dr, M_, options_, udr); end, tElapsed3 = toc(tStart3);

if max(max(abs(ynext1-ynext2))) > 1e-14
    error('Inconsistency between local_state_space_iteration_2 and local_state_space_iteration_k')
end
if max(max(abs(ynext1-ynext3))) > 1e-14
    error('Inconsistency between local_state_space_iteration_2 and local_state_space_iteration_fortran')
end

if tElapsed1<tElapsed2
    skipline()
    dprintf('local_state_space_iteration_2 is %5.2f times faster than local_state_space_iteration_k', tElapsed2/tElapsed1)
    skipline()
else
    skipline()
    dprintf('local_state_space_iteration_2 is %5.2f times slower than local_state_space_iteration_k', tElapsed1/tElapsed2)
    skipline()
end

if tElapsed1<tElapsed3
    skipline()
    dprintf('local_state_space_iteration_2 is %5.2f times faster than local_state_space_iteration_fortran', tElapsed3/tElapsed1)
    skipline()
else
    skipline()
    dprintf('local_state_space_iteration_2 is %5.2f times slower than local_state_space_iteration_fortran', tElapsed1/tElapsed3)
    skipline()
end

stoch_simul(order=3, periods=200, irf=0, k_order_solver);

tStart31 = tic; for i=1:10000, ynext31 = local_state_space_iteration_k(yhat, epsilon, oo_.dr, M_, options_); end, tElapsed31 = toc(tStart31);

tStart32 = tic; [udr] = folded_to_unfolded_dr(oo_.dr, M_, options_); for i=1:10000, ynext32 = local_state_space_iteration_fortran(yhat, epsilon, oo_.dr, M_, options_, udr); end, tElapsed32 = toc(tStart32);

if max(max(abs(ynext31-ynext32))) > 1e-14
    error('Inconsistency between local_state_space_iteration_2 and local_state_space_iteration_fortran')
end

if tElapsed32<tElapsed31
    skipline()
    dprintf('local_state_space_iteration_fortran is %5.2f times faster than local_state_space_iteration_k', tElapsed31/tElapsed32)
    skipline()
else
    skipline()
    dprintf('local_state_space_iteration_fortran is %5.2f times slower than local_state_space_iteration_k', tElapsed32/tElapsed31)
    skipline()
end
