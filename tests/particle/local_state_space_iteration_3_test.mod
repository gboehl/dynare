/*
  Tests that local_state_space_iteration_3 and local_state_space_iteration_k (for k=3) return the same results

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
estimation(datafile='my_data.mat',order=3,mode_compute=0,mh_replic=0,filter_algorithm=sis,nonlinear_filter_initialization=2
    ,cova_compute=0 %tell program that no covariance matrix was computed
);

stoch_simul(order=3, periods=200, irf=0, k_order_solver);

// Really perform the test

nparticles = options_.particle.number_of_particles;
nsims = 1e6/nparticles;

/* We generate particles using realistic distributions (though this is not
   strictly needed) */
nstates = M_.npred + M_.nboth;
state_idx = oo_.dr.order_var((M_.nstatic+1):(M_.nstatic+nstates));
yhat = chol(oo_.var(state_idx,state_idx))*randn(nstates, nparticles);
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
rf_ghxxx = dr.ghxxx(dr.restrict_var_list, :);
rf_ghuuu = dr.ghuuu(dr.restrict_var_list, :);
rf_ghxxu = dr.ghxxu(dr.restrict_var_list, :);
rf_ghxuu = dr.ghxuu(dr.restrict_var_list, :);
rf_ghxss = dr.ghxss(dr.restrict_var_list, :);
rf_ghuss = dr.ghuss(dr.restrict_var_list, :);

options_.threads.local_state_space_iteration_3 = 12;
options_.threads.local_state_space_iteration_k = 12;

% Without pruning
tStart1 = tic; for i=1:nsims, ynext1 = local_state_space_iteration_3(yhat, epsilon, rf_ghx, rf_ghu, rf_constant, rf_ghxx, rf_ghuu, rf_ghxu, rf_ghxxx, rf_ghuuu, rf_ghxxu, rf_ghxuu, rf_ghxss, rf_ghuss, options_.threads.local_state_space_iteration_3); end, tElapsed1 = toc(tStart1);
tStart2 = tic; [udr] = folded_to_unfolded_dr(dr, M_, options_); for i=1:nsims, ynext2 = local_state_space_iteration_k(yhat, epsilon, dr, M_, options_, udr); end, tElapsed2 = toc(tStart2);

if max(max(abs(ynext1-ynext2))) > 1e-10
    error('Without pruning: Inconsistency between local_state_space_iteration_3 and local_state_space_iteration_k')
end

if tElapsed1<tElapsed2
    skipline()
    dprintf('Without pruning: local_state_space_iteration_3 is %5.2f times faster than local_state_space_iteration_k', tElapsed2/tElapsed1)
    skipline()
else
    skipline()
    dprintf('Without pruning: local_state_space_iteration_3 is %5.2f times slower than local_state_space_iteration_k', tElapsed1/tElapsed2)
    skipline()
end

% With pruning
rf_ss = dr.ys(dr.order_var);
rf_ss = rf_ss(dr.restrict_var_list, :);
yhat_ = chol(oo_.var(state_idx,state_idx))*randn(nstates, nparticles);
yhat__ = zeros(2*nstates,nparticles);
yhat__(1:nstates,:) = yhat_;
yhat__(nstates+1:2*nstates,:) = chol(oo_.var(state_idx,state_idx))*randn(nstates, nparticles);
nstatesandobs = size(rf_ghx,1);

[ynext1,ynext1_lat] = local_state_space_iteration_3(yhat, epsilon, rf_ghx, rf_ghu, rf_constant, rf_ghxx, rf_ghuu, rf_ghxu, rf_ghxxx, rf_ghuuu, rf_ghxxu, rf_ghxuu, rf_ghxss, rf_ghuss, yhat__, rf_ss, options_.threads.local_state_space_iteration_3);
[ynext2,ynext2_lat] = local_state_space_iteration_2(yhat__(nstates+1:2*nstates,:), epsilon, rf_ghx, rf_ghu, rf_constant, rf_ghxx, rf_ghuu, rf_ghxu, yhat_, rf_ss, options_.threads.local_state_space_iteration_2);
if max(max(abs(ynext1_lat(nstatesandobs+1:2*nstatesandobs,:)-ynext2))) > 1e-14
    error('With pruning: inconsistency between local_state_space_iteration_2 and local_state_space_iteration_3 on the 2nd-order pruned variable')
end
if max(max(abs(ynext1_lat(1:nstatesandobs,:)-ynext2_lat))) > 1e-14
    error('With pruning: inconsistency between local_state_space_iteration_2 and local_state_space_iteration_3 on the 1st-order pruned variable')
end