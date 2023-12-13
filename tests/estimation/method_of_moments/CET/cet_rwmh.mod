% -------------------------------------------------------------------------
% Functionality testing of Bayesian IRF matching with
% - Random-Walk Metropolis-Hastings
% - loading mode files
% - whether results are stored and can be accessed in different directories
% - reuse previous MCMC to define covariance of proposal
% -------------------------------------------------------------------------

% Copyright © 2023 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

@#include "cet_model.inc"

options_.prior_interval= 0.95;


%%%%%%%%%%%%%%%%%%%%%%
%% NO INTERFACE YET %%
%%%%%%%%%%%%%%%%%%%%%%
[M_.matched_irfs, M_.matched_irfs_weights] = cet_matched_irfs_no_interface_workaround(M_.endo_names,M_.exo_names);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN 1: FIND POSTERIOR MODE USING MODEFILE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method_of_moments(mom_method = irf_matching
%, add_tiny_number_to_cholesky = 1e-14
, additional_optimizer_steps = [4]
%, aim_solver
%, analytic_jacobian
%, analytic_standard_errors
%, bartlett_kernel_lag
%, bounded_shock_support
%, brooks_gelman_plotrows
, cova_compute=0
%, datafile
, dirname=cet_rwmh_results_1
%, dr
%, dr_cycle_reduction_tol
%, dr_logarithmic_reduction_maxiter
%, dr_logarithmic_reduction_tol
%, drop
%, first_obs
%, geweke_interval
%, graph_format
%, huge_number
, irf_matching_file = cet_irf_matching_file
%, k_order_solver
%, load_mh_file
%, load_results_after_load_mh
%, logdata
%, lyapunov
%, lyapunov_complex_threshold
%, lyapunov_doubling_tol
%, lyapunov_fixed_point_tol
%, mcmc_jumping_covariance
, mh_conf_sig = 0.90
%, mh_drop = 0.5
%, mh_init_scale_factor
%, mh_initialize_from_previous_mcmc
%, mh_initialize_from_previous_mcmc_directory
%, mh_initialize_from_previous_mcmc_prior
%, mh_initialize_from_previous_mcmc_record
%, mh_jscale = 0.5
%, mh_nblocks = 1
%, mh_posterior_mode_estimation
%, mh_recover
, mh_replic=0
%, mh_tune_guess = 0.5
%, mh_tune_jscale = 0.33
%, mom_burnin
%, mom_seed
%, mom_se_tolx
, mode_check
%, mode_check_neighbourhood_size
%, mode_check_number_of_points
, mode_check_symmetric_plots = 0
, mode_compute = 1
, mode_file = cet_original_mode
%, nobs
%, no_posterior_kernel_density
%, nodiagnostic
%, nodisplay
%, nograph
%, noprint
%, optim
%, order
%, penalized_estimator
, plot_priors = 1
%, posterior_max_subsample_draws
%, posterior_sampling_method = 'random_walk_metropolis_hastings'
%, posterior_sampler_options = ('proposal_distribution'
%                              ,'rand_multivariate_normal'
%                              ,'rand_multivariate_student'
%                              ,'student_degrees_of_freedom',4
%                              ,'use_mh_covariance_matrix',1
%                              ,'save_tmp_file',1
%                              , scale_file
%                              )
%, posterior_sampling_method = 'tailored_random_block_metropolis_hastings'
%, posterior_sampler_options = ('proposal_distribution'
%                              ,'rand_multivariate_normal'
%                              ,'rand_multivariate_student'
%                              ,'student_degrees_of_freedom',4
%                              ,'new_block_probability',0.3
%                              ,'mode_compute',1
%                              ,optim
%                              ,'save_tmp_file',1
%                              ,scale_file
%                              )
%, posterior_sampling_method = 'independent_metropolis_hastings'
%, posterior_sampler_options = ('proposal_distribution'
%                              ,'rand_multivariate_normal'
%                              ,'rand_multivariate_student'
%                              ,'student_degrees_of_freedom',4
%                              ,'use_mh_covariance_matrix',1
%                              ,'save_tmp_file',1
%                              ,scale_file
%                              )
%, posterior_sampling_method = 'slice'
%, posterior_sampler_options = ('rotated',1
%                              ,'mode_files'
%                              ,'slice_initialize_with_mode',1
%                              ,'initial_step_size',0.8
%                              ,'use_mh_covariance_matrix',1
%                              ,'save_tmp_file',1
%                              )
%, prefilter
%, prior_trunc
%, pruning
%, qz_criterium
%, qz_zero_threshold
%, raftery_lewis_diagnostics
%, raftery_lewis_qrs
%, relative_irf
%, replic
%, schur_vec_tol
%, silent_optimizer
%, simulation_method = stoch_simul
%, simulation_multiple
%, sub_draws
%, sylvester
%, sylvester_fixed_point_tol
%, taper_steps
, tex
%, use_penalized_objective_for_hessian
, verbose
%, weighting_matrix
%, weighting_matrix_scaling_factor
%, xls_range
%, xls_sheet
);
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN 2: LOAD POSTERIOR MODE, COMPUTE HESSIAN AND RUN 1 MCMC CHAIN %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method_of_moments(mom_method = irf_matching
, cova_compute=1
, dirname=cet_rwmh_results_2
, irf_matching_file = cet_irf_matching_file
, mh_conf_sig = 0.95
, mh_drop = 0.5
, mh_jscale = 0.5
, mh_nblocks = 1
, mh_replic=2000
, mode_compute = 0
, mode_file = 'cet_rwmh_results_1/method_of_moments/cet_rwmh_mode'
, plot_priors = 0
, posterior_sampling_method = 'random_walk_metropolis_hastings'
, posterior_sampler_options = ('proposal_distribution'
                              ,'rand_multivariate_normal'
%                              ,'rand_multivariate_student'
%                              ,'student_degrees_of_freedom',4
%                              ,'use_mh_covariance_matrix',1
%                              ,'save_tmp_file',1
%                              , scale_file
                              )
, tex
);
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN 3: DEFINE COVARIANCE OF PROPOSAL FROM PREVIOUS MCMC %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method_of_moments(mom_method = irf_matching
, dirname=cet_rwmh_results_2
, irf_matching_file = cet_irf_matching_file
, load_mh_file
, mh_conf_sig = 0.95
, mh_drop = 0.5
, mh_jscale = 0.5
, mh_nblocks = 1
, mh_replic=600
, mode_compute = 0
, plot_priors = 0
, posterior_sampling_method = 'random_walk_metropolis_hastings'
, posterior_sampler_options = ('proposal_distribution'
%                              ,'rand_multivariate_normal'
                              ,'rand_multivariate_student'
                              ,'student_degrees_of_freedom',4
                              ,'use_mh_covariance_matrix',1
%                              ,'save_tmp_file',1
%                              , scale_file
                              )
, tex
);


%%%%%%%%%%%
%% LATEX %%
%%%%%%%%%%%
collect_latex_files;
if system(['pdflatex -halt-on-error -interaction=batchmode ' M_.fname '_TeX_binder.tex'])
    error('TeX-File did not compile.')
end