@#include "first_spec_common.inc"

varobs q ca;

shocks;
var eeps = 0.04^2;
var nnu = 0.03^2;
var q = 0.01^2;
var ca = 0.01^2;
end;

// Initialize various structures
/*
estimation(datafile='my_data.mat',order=2,mode_compute=0,mh_replic=0,filter_algorithm=sis,nonlinear_filter_initialization=2
    ,cova_compute=0 %tell program that no covariance matrix was computed
);
*/

stoch_simul(order=2, periods=200, irf=0, k_order_solver);

nparticles = 100003;
nsims = 500;
options_.threads.local_state_space_iteration_2 = 32;

state_idx = oo_.dr.order_var((M_.nstatic+1):(M_.nstatic+M_.npred+M_.nboth));

dr = oo_.dr;
dr.restrict_var_list=1:6;


// “rf” stands for “Reduced Form”
rf_ghx = dr.ghx(dr.restrict_var_list, :);
rf_ghu = dr.ghu(dr.restrict_var_list, :);
rf_constant = dr.ys(dr.order_var)+0.5*dr.ghs2;
rf_constant = rf_constant(dr.restrict_var_list, :);
rf_ghxx = dr.ghxx(dr.restrict_var_list, :);
rf_ghuu = dr.ghuu(dr.restrict_var_list, :);
rf_ghxu = dr.ghxu(dr.restrict_var_list, :);

// Check consistency between the kernels
/* We generate particles using realistic distributions (though this is not
   strictly needed) */
yhat = chol(oo_.var(state_idx,state_idx))*randn(M_.npred+M_.nboth, nparticles);
epsilon = chol(M_.Sigma_e)*randn(M_.exo_nbr, nparticles);

setenv('DYNARE_LSSI2_KERNEL', 'avx512')
ynext1 = local_state_space_iteration_2(yhat, epsilon, rf_ghx, rf_ghu, rf_constant, rf_ghxx, rf_ghuu, rf_ghxu, options_.threads.local_state_space_iteration_2);
setenv('DYNARE_LSSI2_KERNEL', 'avx2')
ynext2 = local_state_space_iteration_2(yhat, epsilon, rf_ghx, rf_ghu, rf_constant, rf_ghxx, rf_ghuu, rf_ghxu, options_.threads.local_state_space_iteration_2);
setenv('DYNARE_LSSI2_KERNEL', 'generic')
ynext3 = local_state_space_iteration_2(yhat, epsilon, rf_ghx, rf_ghu, rf_constant, rf_ghxx, rf_ghuu, rf_ghxu, options_.threads.local_state_space_iteration_2);

if max(max(abs(ynext1-ynext3))) > 1e-15
    error('avx512 kernel is inconsistent with generic one')
end
if max(max(abs(ynext2-ynext3))) > 1e-15
    error('avx2 kernel is inconsistent with generic one')
end

// Perform multi-threaded benchmarks
// We regenerate particles before each run to make sure they are not in the CPU cache before the run

tElapsed = NaN(3,nsims);

disp('* Precompute states and shocks for all particles')
yhat_all = NaN(M_.npred+M_.nboth, nparticles, nsims);
epsilon_all = NaN(M_.exo_nbr, nparticles, nsims);
for i = 1:nsims
  yhat_all(:,:,i) = chol(oo_.var(state_idx,state_idx))*randn(M_.npred+M_.nboth, nparticles);
  epsilon_all(:,:,i) = chol(M_.Sigma_e)*randn(M_.exo_nbr, nparticles);
end

disp('* Cool down processor')
pause(10)

disp('* Multi-threaded benchmarks')

for i = 1:nsims

yhat = yhat_all(:,:,i);
epsilon = epsilon_all(:,:,i);

setenv('DYNARE_LSSI2_KERNEL', 'avx512')
tStart = tic;
ynext1 = local_state_space_iteration_2(yhat, epsilon, rf_ghx, rf_ghu, rf_constant, rf_ghxx, rf_ghuu, rf_ghxu, options_.threads.local_state_space_iteration_2);
tElapsed(1,i) = toc(tStart);

setenv('DYNARE_LSSI2_KERNEL', 'avx2')
tStart = tic;
ynext2 = local_state_space_iteration_2(yhat, epsilon, rf_ghx, rf_ghu, rf_constant, rf_ghxx, rf_ghuu, rf_ghxu, options_.threads.local_state_space_iteration_2);
tElapsed(2,i) = toc(tStart);

setenv('DYNARE_LSSI2_KERNEL', 'generic')
tStart = tic;
ynext3 = local_state_space_iteration_2(yhat, epsilon, rf_ghx, rf_ghu, rf_constant, rf_ghxx, rf_ghuu, rf_ghxu, options_.threads.local_state_space_iteration_2);
tElapsed(3,i) = toc(tStart);

end

printf('avx512 kernel:  average timing = %g μs\n', round(mean(tElapsed(1,:))*1e6))
printf('avx2 kernel:    average timing = %g μs\n', round(mean(tElapsed(2,:))*1e6))
printf('generic kernel: average timing = %g μs\n', round(mean(tElapsed(3,:))*1e6))


disp('* Single-threaded benchmarks')

for i = 1:nsims

yhat = yhat_all(:,:,i);
epsilon = epsilon_all(:,:,i);

setenv('DYNARE_LSSI2_KERNEL', 'avx512')
tStart = tic;
ynext1 = local_state_space_iteration_2(yhat, epsilon, rf_ghx, rf_ghu, rf_constant, rf_ghxx, rf_ghuu, rf_ghxu, 1);
tElapsed(1,i) = toc(tStart);

setenv('DYNARE_LSSI2_KERNEL', 'avx2')
tStart = tic;
ynext2 = local_state_space_iteration_2(yhat, epsilon, rf_ghx, rf_ghu, rf_constant, rf_ghxx, rf_ghuu, rf_ghxu, 1);
tElapsed(2,i) = toc(tStart);

setenv('DYNARE_LSSI2_KERNEL', 'generic')
tStart = tic;
ynext3 = local_state_space_iteration_2(yhat, epsilon, rf_ghx, rf_ghu, rf_constant, rf_ghxx, rf_ghuu, rf_ghxu, 1);
tElapsed(3,i) = toc(tStart);

end

printf('avx512 kernel:  average timing = %g μs\n', round(mean(tElapsed(1,:))*1e6))
printf('avx2 kernel:    average timing = %g μs\n', round(mean(tElapsed(2,:))*1e6))
printf('generic kernel: average timing = %g μs\n', round(mean(tElapsed(3,:))*1e6))

setenv('DYNARE_LSSI2_KERNEL', 'auto')
