% -------------------------------------------------------------------------
% Functionality testing of Frequentist IRF matching
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

@#define ML = 1
@#include "cet_model.inc"

options_.prior_interval= 0.95;

method_of_moments(mom_method = irf_matching
%, add_tiny_number_to_cholesky = 1e-14
%, additional_optimizer_steps = [4]
%, aim_solver
%, analytic_jacobian
%, analytic_standard_errors
%, bartlett_kernel_lag
%, bounded_shock_support
%, brooks_gelman_plotrows
%, cova_compute=0
%, datafile
%, dirname=cet_imh_results
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
%, mcmc_jumping_covariance = prior_variance
%, mh_conf_sig = 0.90
%, mh_drop = 0.5
%, mh_init_scale_factor
%, mh_initialize_from_previous_mcmc
%, mh_initialize_from_previous_mcmc_directory
%, mh_initialize_from_previous_mcmc_prior
%, mh_initialize_from_previous_mcmc_record
%, mh_jscale = 0.5
%, mh_nblocks = 2
%, mh_posterior_mode_estimation
%, mh_recover
%, mh_replic=1000
%, mh_tune_guess = 0.5
%, mh_tune_jscale = 0.33
%, mom_burnin
%, mom_seed
%, mom_se_tolx
%, mode_check
%, mode_check_neighbourhood_size
%, mode_check_number_of_points
%, mode_check_symmetric_plots = 0
, mode_compute = 4
%, mode_file = cet_original_mode
%, nobs
%, no_posterior_kernel_density
%, nodiagnostic
%, nodisplay
%, nograph
%, noprint
%, optim
%, order
%, penalized_estimator
%, plot_priors = 0
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
, simulation_method = stoch_simul
%, simulation_multiple
%, sub_draws
%, sylvester
%, sylvester_fixed_point_tol
%, taper_steps
%, tex
%, use_penalized_objective_for_hessian
%, verbose
%, weighting_matrix
%, weighting_matrix_scaling_factor
%, xls_range
%, xls_sheet
);