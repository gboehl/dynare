function [oo_, options_mom_, M_] = run(bayestopt_, options_, oo_, estim_params_, M_, options_mom_)
% function [oo_, options_mom_, M_] = run(bayestopt_, options_, oo_, estim_params_, M_, options_mom_)
% -------------------------------------------------------------------------
% This function performs a method of moments estimation with the following steps:
%  o Checking if required structures and options exist
%  o Preparing local options_mom_ structure
%  o Checking the options and the compatibility of the settings
%  o Initializations of variables, orderings and state space representation
%  o Checks and transformations for matched_moments structure
%  o Checks and transformations for matched_irfs and matched_irfs_weights structure
%  o Checks and transformations for estimated parameters, priors, and bounds
%  o Checks and transformations for data
%  o Checks for objective function at initial parameters
%  o Mode computation: optimization
%    - GMM/SMM: iterated optimization
%    - IRF_MATCHING: optimization
%  o Bayesian MCMC estimation
%  o Display of results
%    - GMM/SMM: J-Test and fit of moments
%    - IRF_MATCHING: fit of IRFs
%  o Clean up
% -------------------------------------------------------------------------
% Note that we call a "mode" the minimum of the objective function, i.e.
% the parameter vector that minimizes the distance between the moments/IRFs
% computed from the model and the moments/IRFs computed from the data.
% -------------------------------------------------------------------------
% This function is inspired by replication codes accompanied to the following papers:
% GMM/SMM:
%  o Andreasen, Fernández-Villaverde, Rubio-Ramírez (2018): "The Pruned State-Space System for Non-Linear DSGE Models: Theory and Empirical Applications", Review of Economic Studies, 85(1):1-49.
%  o Born, Pfeifer (2014): "Risk Matters: Comment", American Economic Review, 104(12):4231-4239.
%  o Mutschler (2018): "Higher-order statistics for DSGE models", Econometrics and Statistics, 6:44-56.
%  o Ruge-Murcia (2007): "Methods to Estimate Dynamic Stochastic General Equilibrium Models", Journal of Economic Dynamics and Control, 31(8):2599-2636.
% IRF MATCHING:
%  o Christiano, Trabandt, Walentin (2010): "DSGE Models for Monetary Policy Analysis." In Handbook of Monetary Economics, 3:285–367.
%  o Christiano, Eichenbaum, Trabandt (2016): "Unemployment and Business Cycles." Econometrica, 84: 1523-1569.
%  o Ruge-Murcia (2020): "Estimating Nonlinear Dynamic Equilibrium Models by Matching Impulse Responses", Economics Letters, 197.
% -------------------------------------------------------------------------
% INPUTS
%  o bayestopt_:     [structure] information about priors
%  o options_:       [structure] information about global options
%  o oo_:            [structure] results
%  o estim_params_:  [structure] information about estimated parameters
%  o M_:             [structure] information about model with
%                     o matched_moments:  [cell] information about selected moments to match in GMM/SMM estimation
%                                                vars: matched_moments{:,1});
%                                                lead/lags: matched_moments{:,2};
%                                                powers: matched_moments{:,3};
%                     o matched_irfs:          [cell] information about selected IRFs to match in IRF_MATCHING estimation
%                     o matched_irfs_weights:  [cell] information about entries in weight matrix for an IRF_MATCHING estimation
%  o options_mom_:   [structure] information about settings specified by the user
% -------------------------------------------------------------------------
% OUTPUTS
%  o oo_:            [structure] storage for results (oo_)
%  o options_mom_:   [structure] information about all (user-specified and updated) settings used in estimation (options_mom_)
%  o M_:             [structure] updated information about model
% -------------------------------------------------------------------------
% This function is called by
%  o driver.m
% -------------------------------------------------------------------------
% This function calls
%  o check_for_calibrated_covariances
%  o check_mode_file
%  o check_posterior_sampler_options
%  o check_prior_bounds
%  o check_prior_stderr_corr
%  o check_steady_state_changes_parameters
%  o check_varobs_are_endo_and_declared_once
%  o check_hessian_at_the_mode
%  o display_estimation_results_table
%  o do_parameter_initialization
%  o get_all_parameters
%  o get_dynare_random_generator_state
%  o get_matrix_entries_for_psd_check
%  o M_.fname '_prior_restrictions'
%  o makedataset
%  o mode_check
%  o mom.check_irf_matching_file
%  o mom.default_option_mom_values
%  o mom.get_data_moments
%  o mom.matched_irfs_blocks
%  o mom.matched_moments_block
%  o mom.objective_function
%  o mom.optimal_weighting_matrix
%  o mom.print_info_on_estimation_settings
%  o mom.set_correct_bounds_for_stderr_corr
%  o mom.standard_errors
%  o plot_priors
%  o prior_bounds
%  o priordens
%  o print_info
%  o set_all_parameters
%  o set_dynare_random_generator_state
%  o set_prior
%  o set_state_space
%  o skipline
%  o test_for_deep_parameters_calibration
%  o transform_prior_to_laplace_prior

% Copyright © 2020-2023 Dynare Team
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

fprintf('\n==== Method of Moments Estimation (%s) ====\n\n',options_mom_.mom.mom_method);


% -------------------------------------------------------------------------
% checks if required structures exist
% -------------------------------------------------------------------------
if isempty(estim_params_) % structure storing the info about estimated parameters in the estimated_params block
    if ~(isfield(estim_params_,'nvx') && (size(estim_params_.var_exo,1)+size(estim_params_.var_endo,1)+size(estim_params_.corrx,1)+size(estim_params_.corrn,1)+size(estim_params_.param_vals,1))==0)
        error('method_of_moments: You need to provide an ''estimated_params'' block!');
    else
        error('method_of_moments: The ''estimated_params'' block must not be empty!');
    end
end
if strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')
    if ~isfield(M_,'matched_moments') || isempty(M_.matched_moments) % structure storing the moments used for GMM and SMM estimation
        error('method_of_moments: You need to provide a ''matched_moments'' block for ''mom_method=%s''!',options_mom_.mom.mom_method);
    end
elseif strcmp(options_mom_.mom.mom_method,'IRF_MATCHING')
    if ~isfield(M_,'matched_irfs') || isempty(M_.matched_irfs) % structure storing the irfs used for matching
        error('method_of_moments: You need to provide a ''matched_irfs'' block for ''mom_method=%s''!',options_mom_.mom.mom_method);
    end
end
if (~isempty(estim_params_.var_endo) || ~isempty(estim_params_.corrn)) && strcmp(options_mom_.mom.mom_method, 'GMM')
    error('method_of_moments: GMM estimation does not support measurement error(s) yet. Please specify them as a structural shock!');
end
do_bayesian_estimation = [estim_params_.var_exo(:,5); estim_params_.var_endo(:,5); estim_params_.corrx(:,6); estim_params_.corrn(:,6); estim_params_.param_vals(:,5)];
if all(do_bayesian_estimation~=0)
    do_bayesian_estimation = true;
elseif all(do_bayesian_estimation==0)
    do_bayesian_estimation = false;
else
    error('method_of_moments: Estimation must be either fully Frequentist or fully Bayesian. Maybe you forgot to specify a prior distribution!');
end
if ~isfield(options_,'varobs')
    error('method_of_moments: VAROBS statement is missing!');
end
check_varobs_are_endo_and_declared_once(options_.varobs,M_.endo_names);


% -------------------------------------------------------------------------
% options_mom_ structure
% -------------------------------------------------------------------------
% options_mom_ is local and contains default and user-specified values for
% all settings needed for the method of moments estimation. Some options,
% though, are set by the preprocessor into options_ and we copy these over.
% The idea is to be independent of options_ and have full control of the
% estimation instead of possibly having to deal with options chosen somewhere
% else in the mod file.
options_mom_ = mom.default_option_mom_values(options_mom_, options_, M_.dname, M_.fname, do_bayesian_estimation);


% -------------------------------------------------------------------------
% workarounds
% -------------------------------------------------------------------------
if strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')
% temporary workaround for https://git.dynare.org/Dynare/dseries/-/issues/51
    if options_mom_.xls_sheet~=1
        evalin('base','options_.xls_sheet=options_mom_.xls_sheet');
    end
    if ~isempty(options_mom_.xls_range)
        evalin('base','options_.xls_range=options_mom_.xls_range');
    end
end


% -------------------------------------------------------------------------
% initializations
% -------------------------------------------------------------------------
%save warning state for restoring later on
orig_warning_state = warning;
% create output directories to store results
M_.dname = options_mom_.dirname;
CheckPath(M_.dname,'.');
CheckPath('method_of_moments',M_.dname);
CheckPath('graphs',M_.dname);

if do_bayesian_estimation
    oo_.mom.posterior.optimization.mode = [];
    oo_.mom.posterior.optimization.Variance = [];
    oo_.mom.posterior.optimization.log_density=[];
end
do_bayesian_estimation_mcmc = do_bayesian_estimation && ( (options_mom_.mh_replic>0) || options_mom_.load_mh_file );
invhess = [];
% decision rule
oo_.dr = set_state_space(oo_.dr,M_); % get state-space representation
options_mom_.mom.obs_var = []; % create index of observed variables in DR order
for i = 1:options_mom_.obs_nbr
    options_mom_.mom.obs_var = [options_mom_.mom.obs_var; find(strcmp(options_mom_.varobs{i}, M_.endo_names(oo_.dr.order_var)))];
end


% -------------------------------------------------------------------------
% matched_moments: checks and transformations
% -------------------------------------------------------------------------
if strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')
    M_.matched_moments = mom.matched_moments_block(M_.matched_moments, options_mom_.mom.mom_method);    
    % Check if both prefilter and first moments were specified
    first_moment_indicator = find(cellfun(@(x) sum(abs(x))==1,M_.matched_moments(:,3)))';
    if options_mom_.prefilter && ~isempty(first_moment_indicator)
        fprintf('Centered moments requested (prefilter option is set); therefore, ignore declared first moments in ''matched_moments'' block.\n');
        M_.matched_moments(first_moment_indicator,:)=[]; %remove first moments entries
    end
    options_mom_.mom.mom_nbr = size(M_.matched_moments,1);
    % Get maximum lag number for autocovariances/autocorrelations
    options_mom_.ar = max(cellfun(@max,M_.matched_moments(:,2))) - min(cellfun(@min,M_.matched_moments(:,2)));
    % Check that only observed variables are involved in moments
    not_observed_variables=setdiff(oo_.dr.inv_order_var([M_.matched_moments{:,1}]),options_mom_.mom.obs_var);
    if ~isempty(not_observed_variables)
        skipline;
        error('method_of_moments: You specified moments involving %s, but it is not a varobs!',M_.endo_names{oo_.dr.order_var(not_observed_variables)});
    end
end


% -------------------------------------------------------------------------
% matched_irfs: checks and transformations
% -------------------------------------------------------------------------
if strcmp(options_mom_.mom.mom_method,'IRF_MATCHING')
    [oo_.mom.data_moments, oo_.mom.weighting_info.W, options_mom_.mom.irfIndex, options_mom_.irf] = mom.matched_irfs_blocks(M_.matched_irfs, M_.matched_irfs_weights, options_mom_.varobs_id, options_mom_.obs_nbr, M_.exo_nbr, M_.endo_names, M_.exo_names);
    % compute inverse of weighting matrix
    try
        oo_.mom.weighting_info.Winv = inv(oo_.mom.weighting_info.W);
    catch
        error('method_of_moments: Something wrong while computing inv(W), check your weighting matrix!');
    end
    if any(isnan(oo_.mom.weighting_info.Winv(:))) || any(isinf(oo_.mom.weighting_info.Winv(:)))
        error('method_of_moments: There are NaN or Inf values in inv(W), check your weighting matrix!');
    end
    % compute log determinant of inverse of weighting matrix in a robust way to avoid Inf or NaN
    try
        oo_.mom.weighting_info.Winv_logdet = 2*sum(log(diag(chol(oo_.mom.weighting_info.Winv))));
    catch
        error('method_of_moments: Something wrong while computing log(det(inv(W))), check your weighting matrix!');
    end
    if any(isnan(oo_.mom.weighting_info.Winv_logdet(:))) || any(isinf(oo_.mom.weighting_info.Winv_logdet(:)))
        error('method_of_moments: There are NaN or Inf values in log(det(inv(W))), check your weighting matrix!');
    end
    options_mom_.mom.mom_nbr = length(options_mom_.mom.irfIndex);
end


% -------------------------------------------------------------------------
% irf_matching_file: checks and transformations
% -------------------------------------------------------------------------
if strcmp(options_mom_.mom.mom_method,'IRF_MATCHING')
    [options_mom_.mom.irf_matching_file.name, options_mom_.mom.irf_matching_file.path] = mom.check_irf_matching_file(options_mom_.mom.irf_matching_file.name);
    % check for irf_matching_file
    if ~( isempty(options_mom_.mom.irf_matching_file.path) || strcmp(options_mom_.mom.irf_matching_file.path,'.') )
        fprintf('\nAdding %s to MATLAB''s path.\n',options_mom_.mom.irf_matching_file.path);
        addpath(options_mom_.mom.irf_matching_file.path);
    end
end


% -------------------------------------------------------------------------
% estimated parameters: checks and transformations on values, priors, bounds, posterior options
% -------------------------------------------------------------------------
% set priors and bounds over the estimated parameters
[xparam0, estim_params_, bayestopt_, lb, ub, M_] = set_prior(estim_params_, M_, options_mom_);
number_of_estimated_parameters = length(xparam0);
hessian_xparam0 = []; % initialize hessian
% check if enough moments for estimation
if options_mom_.mom.mom_nbr < length(xparam0)
    skipline;
    error('method_of_moments: There must be at least as many moments as parameters for a %s estimation!',options_mom_.mom.mom_method);
end
skipline(2);
% check if a _prior_restrictions.m file exists
if exist([M_.fname '_prior_restrictions.m'],'file')
    options_mom_.prior_restrictions.status = 1;
    options_mom_.prior_restrictions.routine = str2func([M_.fname '_prior_restrictions']);
end
% check that the provided mode_file is compatible with the current estimation settings
if ~isempty(options_mom_.mode_file) && ( ~do_bayesian_estimation || (do_bayesian_estimation && ~options_mom_.mh_posterior_mode_estimation) )
    [xparam0, hessian_xparam0] = check_mode_file(xparam0, hessian_xparam0, options_mom_, bayestopt_);
end
% check on specified priors and penalized estimation (which uses Laplace approximated priors)
if strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')
    bayestopt_orig = bayestopt_;
    if any(bayestopt_.pshape > 0) % prior specified
        if ~options_mom_.mom.penalized_estimator
            fprintf('\nPriors were specified, but the penalized_estimator-option was not set.');
            fprintf('\nDynare sets penalized_estimator to 1. Conducting %s with penalty.\n',options_mom_.mom.mom_method);
            options_mom_.mom.penalized_estimator = 1;
        end
        bayestopt_ = mom.transform_prior_to_laplace_prior(bayestopt_);
    end
end
% check for calibrated covariances before updating parameters
estim_params_ = check_for_calibrated_covariances(estim_params_,M_);

% checks on parameter calibration and initialization
xparam_calib = get_all_parameters(estim_params_,M_); % get calibrated parameters
if ~any(isnan(xparam_calib)) % all estimated parameters are calibrated
    estim_params_.full_calibration_detected = true;
else
    estim_params_.full_calibration_detected = false;
end
if options_mom_.use_calibration_initialization % set calibration as starting values
    if ~isempty(bayestopt_) && ~do_bayesian_estimation && any(all(isnan([xparam_calib xparam0]),2))
        error('method_of_moments: When using the use_calibration option with %s without prior, the parameters must be explicitly initialized!',options_mom_.mom.mom_method);
    else
        [xparam0,estim_params_] = do_parameter_initialization(estim_params_,xparam_calib,xparam0); % get explicitly initialized parameters that have precedence over calibrated values
    end
end
% check initialization
if ~isempty(bayestopt_) && ~do_bayesian_estimation && any(isnan(xparam0))
    error('method_of_moments: Frequentist %s requires all estimated parameters to be initialized, either in an estimated_params or estimated_params_init-block!',options_mom_.mom.mom_method);
end
% set and check parameter bounds
if ~isempty(bayestopt_) && do_bayesian_estimation
    % plot prior densities
    if ~options_mom_.nograph && options_mom_.plot_priors
        if strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')
            plot_priors(bayestopt_orig,M_,estim_params_,options_mom_,'Original priors'); % only for visual inspection (not saved to disk, because overwritten in next call to plot_priors)
            plot_priors(bayestopt_,M_,estim_params_,options_mom_,'Laplace approximated priors');
            clear('bayestopt_orig'); % make sure stale structure cannot be used
        else
            plot_priors(bayestopt_,M_,estim_params_,options_mom_,'Priors');
        end
    end
    % set prior bounds
    BoundsInfo = prior_bounds(bayestopt_, options_mom_.prior_trunc);
    BoundsInfo.lb = max(BoundsInfo.lb,lb);
    BoundsInfo.ub = min(BoundsInfo.ub,ub);
else
    % no priors are declared so Dynare will estimate the parameters with Frequentist methods using inequality constraints for the parameters
    BoundsInfo.lb = lb;
    BoundsInfo.ub = ub;
    if (strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')) && options_mom_.mom.penalized_estimator
        fprintf('Penalized estimation turned off as you did not declare priors\n');
        options_mom_.mom.penalized_estimator = 0;
    else
        if isfield(options_mom_,'mh_replic') && options_mom_.mh_replic > 0
            fprintf('Setting ''mh_replic=0'' as you did not declare priors.\n');
            options_mom_.mh_replic = 0;
        end
    end
end

% set correct bounds for standard deviations and correlations
BoundsInfo = mom.set_correct_bounds_for_stderr_corr(estim_params_,BoundsInfo);
% test if initial values of the estimated parameters are all between the prior lower and upper bounds
if options_mom_.use_calibration_initialization
    try
        check_prior_bounds(xparam0,BoundsInfo,M_,estim_params_,options_mom_,bayestopt_);
    catch last_error
        fprintf('Cannot use parameter values from calibration as they violate the prior bounds.');
        rethrow(last_error);
    end
else
    check_prior_bounds(xparam0,BoundsInfo,M_,estim_params_,options_mom_,bayestopt_);
end
% check for positive definiteness
estim_params_ = get_matrix_entries_for_psd_check(M_,estim_params_);
% set sigma_e_is_diagonal flag (needed if the shocks block is not declared in the mod file)
M_.sigma_e_is_diagonal = true;
if estim_params_.ncx || any(nnz(tril(M_.Correlation_matrix,-1))) || isfield(estim_params_,'calibrated_covariances')
    M_.sigma_e_is_diagonal = false;
end
% storing prior parameters in results structure
if do_bayesian_estimation || ( (strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')) && options_mom_.mom.penalized_estimator)
    oo_.mom.prior.mean = bayestopt_.p1;
    oo_.mom.prior.mode = bayestopt_.p5;
    oo_.mom.prior.variance = diag(bayestopt_.p2.^2);
    oo_.mom.prior.hyperparameters.first = bayestopt_.p6;
    oo_.mom.prior.hyperparameters.second = bayestopt_.p7;
end
% set all parameters
M_ = set_all_parameters(xparam0,estim_params_,M_);
% provide warning if there is NaN in parameters
test_for_deep_parameters_calibration(M_);
% set jscale
if do_bayesian_estimation_mcmc
    if ~strcmp(options_mom_.posterior_sampler_options.posterior_sampling_method,'slice')
        if isempty(options_mom_.mh_jscale)
            options_mom_.mh_jscale = 2.38/sqrt(number_of_estimated_parameters); % use optimal value for univariate normal distribution (check_posterior_sampler_options and mode_compute=6 may overwrite this setting)
        end
        bayestopt_.jscale(find(isnan(bayestopt_.jscale))) = options_mom_.mh_jscale;
    end
end
% initialization of posterior sampler options
if do_bayesian_estimation_mcmc
    [current_options, options_mom_, bayestopt_] = check_posterior_sampler_options([], M_.fname, M_.dname, options_mom_, BoundsInfo, bayestopt_);
    options_mom_.posterior_sampler_options.current_options = current_options;
    if strcmp(current_options.posterior_sampling_method,'slice') && current_options.use_mh_covariance_matrix && ~current_options.rotated
        error('method_of_moments: Using the slice sampler with the ''use_mh_covariance_matrix'' option requires also setting the ''rotated'' option!');
    end
end
% warning if prior allows that stderr parameters are negative or corr parameters are outside the unit circle
if do_bayesian_estimation
    check_prior_stderr_corr(estim_params_,bayestopt_);
    % check value of prior density
    [~,~,~,info] = priordens(xparam0,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4);
    if any(info)
        fprintf('The prior density evaluated at the initial values is Inf for the following parameters: %s\n',bayestopt_.name{info,1});
        error('The initial value of the prior is -Inf!');
    end
end


% -------------------------------------------------------------------------
% datafile: checks and transformations
% -------------------------------------------------------------------------
% build dataset
if strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')
    % check if datafile has same name as mod file
    [~,name] = fileparts(options_mom_.datafile);
    if strcmp(name,M_.fname)
        error('method_of_moments: ''datafile'' and mod file are not allowed to have the same name; change the name of the ''datafile''!');
    end
    dataset_ = makedataset(options_mom_);
    % set options for old interface from the ones for new interface
    if ~isempty(dataset_)
        options_mom_.nobs = dataset_.nobs;
    end
    % check length of data for estimation of second moments
    if options_mom_.ar > options_mom_.nobs+1
        error('method_of_moments: Dataset is too short to compute higher than first moments!');
    end
    % provide info on data moments handling
    fprintf('Computing data moments. Note that NaN values in the moments (due to leads and lags or missing data) are replaced by the mean of the corresponding moment.\n');
    % get data moments for the method of moments
    [oo_.mom.data_moments, oo_.mom.m_data] = mom.get_data_moments(dataset_.data, options_mom_.mom.obs_var, oo_.dr.inv_order_var, M_.matched_moments, options_mom_);
    if ~isreal(dataset_.data)
        error('method_of_moments: The data moments contain complex values!');
    end
end


% -------------------------------------------------------------------------
% SMM: Get shock series and set variance correction factor
% -------------------------------------------------------------------------
if strcmp(options_mom_.mom.mom_method,'SMM')
    options_mom_.mom.long = round(options_mom_.mom.simulation_multiple*options_mom_.nobs);
    options_mom_.mom.variance_correction_factor = (1+1/options_mom_.mom.simulation_multiple);
    % draw shocks for SMM
    if ~isoctave
        smmstream = RandStream('mt19937ar','Seed',options_mom_.mom.seed);
        temp_shocks = randn(smmstream,options_mom_.mom.long+options_mom_.mom.burnin,M_.exo_nbr);
        temp_shocks_ME = randn(smmstream,options_mom_.mom.long,length(M_.H));
    else
        [state_u,state_n] = get_dynare_random_generator_state; %get state for later resetting
        set_dynare_random_generator_state(options_mom_.mom.seed,options_mom_.mom.seed);
        temp_shocks = randn(options_mom_.mom.long+options_mom_.mom.burnin,M_.exo_nbr);
        temp_shocks_ME = randn(options_mom_.mom.long,length(M_.H));
        set_dynare_random_generator_state(state_u,state_n); %reset state for later resetting
    end
    if options_mom_.mom.bounded_shock_support == 1
        temp_shocks(temp_shocks>2) = 2;
        temp_shocks(temp_shocks<-2) = -2;
        temp_shocks_ME(temp_shocks_ME<-2) = -2;
        temp_shocks_ME(temp_shocks_ME<-2) = -2;
    end
    options_mom_.mom.shock_series = temp_shocks;
    options_mom_.mom.ME_shock_series = temp_shocks_ME;
    if options_mom_.k_order_solver && ~options_mom_.pruning % dynare++ routines will be called in simult_.m, store some additional stuff
        options_mom_.DynareRandomStreams.seed = options_mom_.mom.seed;
    end
end


% -------------------------------------------------------------------------
% checks for steady state at initial parameters
% -------------------------------------------------------------------------
% check if steady state solves static model and if steady-state changes estimated parameters
if options_mom_.steadystate.nocheck
    steadystate_check_flag_vec = [0 1];
else
    steadystate_check_flag_vec = [1 1];
end
[oo_.steady_state, info, steady_state_changes_parameters] = check_steady_state_changes_parameters(M_, estim_params_, oo_, options_mom_, steadystate_check_flag_vec);
if info(1)
    fprintf('\nThe steady state at the initial parameters cannot be computed.\n');
    print_info(info, 0, options_mom_);
end
if steady_state_changes_parameters && strcmp(options_mom_.mom.mom_method,'GMM') && options_mom_.mom.analytic_standard_errors
    fprintf('For analytical standard errors, the parameter-Jacobians of the dynamic model and of the steady-state will be computed numerically,\n');
    fprintf('because the steady-state changes estimated parameters. Option ''analytic_derivation_mode'' reset to -2.');
    options_mom_.analytic_derivation_mode = -2;
end
% display warning if some parameters are still NaN
test_for_deep_parameters_calibration(M_);


% -------------------------------------------------------------------------
% checks for objective function at initial parameters
% -------------------------------------------------------------------------
objective_function = str2func('mom.objective_function');
try
    % check for NaN or complex values of moment-distance-funtion evaluated at initial parameters
    if strcmp(options_mom_.mom.mom_method,'SMM') || strcmp(options_mom_.mom.mom_method,'GMM')
        oo_.mom.weighting_info.Sw = eye(options_mom_.mom.mom_nbr); % initialize with identity weighting matrix
    end
    tic_id = tic;
    [fval, info] = feval(objective_function, xparam0, oo_.mom.data_moments, oo_.mom.weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
    elapsed_time = toc(tic_id);
    if isnan(fval)
        if strcmp(options_mom_.mom.mom_method,'SMM') || strcmp(options_mom_.mom.mom_method,'GMM')
            error('method_of_moments: The initial value of the objective function with identity weighting matrix is NaN!');
        else
            error('method_of_moments: The initial value of the objective function is NaN!');
        end
    elseif imag(fval)
        if strcmp(options_mom_.mom.mom_method,'SMM') || strcmp(options_mom_.mom.mom_method,'GMM')
            error('method_of_moments: The initial value of the objective function with identity weighting matrix is complex!');
        else
            error('method_of_moments: The initial value of the objective function is complex!');
        end
    end
    if info(1) > 0
        disp('method_of_moments: Error in computing the objective function for initial parameter values')
        print_info(info, options_mom_.noprint, options_mom_)
    end
    fprintf('Initial value of the moment objective function');
    if strcmp(options_mom_.mom.mom_method,'SMM') || strcmp(options_mom_.mom.mom_method,'GMM')
        fprintf(' with %4.1f times identity weighting matrix', options_mom_.mom.weighting_matrix_scaling_factor);
    end
    fprintf(': %6.4f \n\n', fval);
    fprintf('Time required to compute objective function once: %5.4f seconds \n', elapsed_time);
catch last_error % if check fails, provide info on using calibration if present
    if estim_params_.full_calibration_detected %calibrated model present and no explicit starting values
        skipline(1);
        fprintf('There was an error in computing the moments for initial parameter values.\n');
        fprintf('If this is not a problem with the setting of options (check the error message below),\n');
        fprintf('you should try using the calibrated version of the model as starting values. To do\n');
        fprintf('this, add an empty estimated_params_init-block with use_calibration option immediately before the estimation\n');
        fprintf('command (and after the estimated_params-block so that it does not get overwritten):\n');
        skipline(2);
    end
    rethrow(last_error);
end


% -------------------------------------------------------------------------
% print some info to console
% -------------------------------------------------------------------------
mom.print_info_on_estimation_settings(options_mom_, number_of_estimated_parameters, do_bayesian_estimation);


% -------------------------------------------------------------------------
% compute mode for GMM/SMM
% -------------------------------------------------------------------------
if strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')
    [xparam1, oo_.mom.weighting_info, oo_.mom.verbose] = mom.mode_compute_gmm_smm(xparam0, objective_function, oo_.mom.m_data, oo_.mom.data_moments, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
end

% -------------------------------------------------------------------------
% compute mode for IRF matching
% -------------------------------------------------------------------------
if strcmp(options_mom_.mom.mom_method,'IRF_MATCHING')
    if ~do_bayesian_estimation || (do_bayesian_estimation && ~options_mom_.mh_posterior_mode_estimation)
        [xparam1, hessian_xparam1, fval, oo_.mom.verbose] = mom.mode_compute_irf_matching(xparam0, hessian_xparam0, objective_function, do_bayesian_estimation, oo_.mom.weighting_info, oo_.mom.data_moments, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
    else
        xparam1 = xparam0;
        hessian_xparam1 = hessian_xparam0;
    end
end


% -------------------------------------------------------------------------
% compute standard errors and initialize covariance of the proposal distribution
% -------------------------------------------------------------------------
if strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')
    % compute standard errors at mode
    options_mom_.mom.vector_output = false; % make sure flag is reset
    M_ = set_all_parameters(xparam1,estim_params_,M_); % update M_ and oo_ (in particular to get oo_.mom.model_moments)
    if strcmp(options_mom_.mom.mom_method,'GMM') && options_mom_.mom.analytic_standard_errors
        options_mom_.mom.compute_derivs = true; % for GMM we compute derivatives analytically in the objective function with this flag
    end
    [~, ~, ~, ~, ~, oo_.mom.Q, oo_.mom.model_moments, oo_.mom.model_moments_params_derivs] = feval(objective_function, xparam1, oo_.mom.data_moments, oo_.mom.weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
    options_mom_.mom.compute_derivs = false; % reset to not compute derivatives in objective function during optimization
    [stdh, invhess] = mom.standard_errors(xparam1, objective_function, oo_.mom.model_moments, oo_.mom.model_moments_params_derivs, oo_.mom.m_data, oo_.mom.data_moments, oo_.mom.weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
    if options_mom_.cova_compute
        hessian_xparam1 = inv(invhess);
    end
else
    if ~do_bayesian_estimation || ~options_mom_.mh_posterior_mode_estimation
        if do_bayesian_estimation
            oo_.mom.posterior.optimization.mode = xparam1;
            if exist('fval','var')
                oo_.mom.posterior.optimization.log_density = -fval;
            end
        end
        if options_mom_.cova_compute
            hsd = sqrt(diag(hessian_xparam1)); % represent the curvature (or second derivatives) of the likelihood with respect to each parameter being estimated.
            invhess = inv(hessian_xparam1./(hsd*hsd'))./(hsd*hsd'); % before taking the inverse scale the Hessian matrix by dividing each of its elements by the outer product of hsd such that the diagonal of the resulting matrix is approximately 1. This kind of scaling can help in regularizing the matrix and potentially improves its condition number, which in turn can make the matrix inversion more stable.
            stdh = sqrt(diag(invhess));
            if do_bayesian_estimation
                oo_.mom.posterior.optimization.Variance = invhess;
            end
        end
    else
        variances = bayestopt_.p2.*bayestopt_.p2;
        id_Inf = isinf(variances);
        variances(id_Inf) = 1;
        invhess = options_mom_.mh_posterior_mode_estimation*diag(variances);
        xparam1 = bayestopt_.p5;
        id_NaN = isnan(xparam1);
        xparam1(id_NaN) = bayestopt_.p1(id_NaN);
        outside_bound_pars=find(xparam1 < BoundsInfo.lb | xparam1 > BoundsInfo.ub);
        xparam1(outside_bound_pars) = bayestopt_.p1(outside_bound_pars);
    end
    if ~options_mom_.cova_compute
        stdh = NaN(length(xparam1),1);
    end
end


% -------------------------------------------------------------------------
% display estimation results at mode
% -------------------------------------------------------------------------
if do_bayesian_estimation && ~options_mom_.mom.penalized_estimator && ~options_mom_.mh_posterior_mode_estimation
    % display table with Bayesian mode estimation results and store parameter estimates and standard errors in oo_
    oo_.mom = display_estimation_results_table(xparam1, stdh, M_, options_mom_, estim_params_, bayestopt_, oo_.mom, prior_dist_names, 'Posterior', 'posterior');
    % Laplace approximation to the marginal log density
    if options_mom_.cova_compute
        estim_params_nbr = size(xparam1,1);
        if ispd(invhess)
            log_det_invhess = log(det(invhess./(stdh*stdh')))+2*sum(log(stdh));
            likelihood = feval(objective_function, xparam1, oo_.mom.data_moments, oo_.mom.weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
            oo_.mom.MarginalDensity.LaplaceApproximation = .5*estim_params_nbr*log(2*pi) + .5*log_det_invhess - likelihood;
        else
            oo_.mom.MarginalDensity.LaplaceApproximation = NaN;
        end        
        fprintf('\nLog data density [Laplace approximation] is %f.\n',oo_.mom.MarginalDensity.LaplaceApproximation);
    end
elseif ~do_bayesian_estimation || (do_bayesian_estimation && options_mom_.mom.penalized_estimator)
    % display table with Frequentist estimation results and store parameter estimates and standard errors in oo_
    oo_.mom = display_estimation_results_table(xparam1, stdh, M_, options_mom_, estim_params_, bayestopt_, oo_.mom, prior_dist_names, options_mom_.mom.mom_method, lower(options_mom_.mom.mom_method));    
end


% -------------------------------------------------------------------------
% checks for mode and hessian at the mode
% -------------------------------------------------------------------------
if (~do_bayesian_estimation && options_mom_.cova_compute) || (do_bayesian_estimation && ~options_mom_.mh_posterior_mode_estimation && options_mom_.cova_compute)
    check_hessian_at_the_mode(hessian_xparam1, xparam1, M_, estim_params_, options_, BoundsInfo);
end
if options_mom_.mode_check.status
    if ~do_bayesian_estimation || (do_bayesian_estimation && ~options_mom_.mh_posterior_mode_estimation)
        mode_check(objective_function, xparam1, diag(stdh), options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, true,... % use diag(stdh) instead of hessian_xparam1 as mode_check uses diagonal elements
                   oo_.mom.data_moments, oo_.mom.weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
    end
end

% -------------------------------------------------------------------------
% Bayesian MCMC estimation
% -------------------------------------------------------------------------
if do_bayesian_estimation_mcmc
    invhess = set_mcmc_jumping_covariance(invhess, length(xparam1), options_mom_.MCMC_jumping_covariance, bayestopt_, 'method_of_moments');
    % reset bounds as lb and ub must only be operational during mode-finding
    BoundsInfo = set_mcmc_prior_bounds(xparam1, bayestopt_, options_mom_, 'method_of_moments');
    % tunes the jumping distribution's scale parameter
    if isfield(options_mom_,'mh_tune_jscale') && options_mom_.mh_tune_jscale.status
        if strcmp(options_mom_.posterior_sampler_options.posterior_sampling_method, 'random_walk_metropolis_hastings')
            options_mom_.mh_jscale = tune_mcmc_mh_jscale_wrapper(invhess, options_mom_, M_, objective_function, xparam1, BoundsInfo,...
                                                                 oo_.mom.data_moments, oo_.mom.weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
            bayestopt_.jscale(:) = options_mom_.mh_jscale;
            fprintf('mh_tune_jscale: mh_jscale has been set equal to %s.\n', num2str(options_mom_.mh_jscale));
        else
            warning('mh_tune_jscale is only available with ''random_walk_metropolis_hastings''!')
        end
    end
    % run MCMC sampling
    posterior_sampler_options = options_mom_.posterior_sampler_options.current_options;
    posterior_sampler_options.invhess = invhess;
    [posterior_sampler_options, options_mom_, bayestopt_] = check_posterior_sampler_options(posterior_sampler_options, M_.fname, M_.dname, options_mom_, BoundsInfo, bayestopt_,'method_of_moments');    
    options_mom_.posterior_sampler_options.current_options = posterior_sampler_options; % store current options
    if options_mom_.mh_replic>0
        posterior_sampler(objective_function,posterior_sampler_options.proposal_distribution,xparam1,posterior_sampler_options,BoundsInfo,oo_.mom.data_moments,oo_.mom.weighting_info,options_mom_,M_,estim_params_,bayestopt_,oo_,'method_of_moments::mcmc');
    end
    CutSample(M_, options_mom_, 'method_of_moments::mcmc'); % discard first mh_drop percent of the draws
    if options_mom_.mh_posterior_mode_estimation
        % skip optimizer-based mode-finding and instead compute the mode based on a run of a MCMC
        [~,~,posterior_mode,~] = compute_mh_covariance_matrix(bayestopt_,M_.fname,M_.dname,'method_of_moments');
        oo_.mom = fill_mh_mode(posterior_mode',NaN(length(posterior_mode),1),M_,options_mom_,estim_params_,bayestopt_,oo_.mom,'posterior');
        warning(orig_warning_state);
        return
    else
        % get stored results if required
        if options_mom_.load_mh_file && options_mom_.load_results_after_load_mh
            oo_load_mh = load([M_.dname filesep 'method_of_moments' filesep M_.fname '_mom_results'],'oo_');
        end
        % convergence diagnostics
        if ~options_mom_.nodiagnostic
            if (options_mom_.mh_replic>0 || (options_mom_.load_mh_file && ~options_mom_.load_results_after_load_mh))
                oo_.mom = mcmc_diagnostics(options_mom_, estim_params_, M_, oo_.mom);
            elseif options_mom_.load_mh_file && options_mom_.load_results_after_load_mh
                if isfield(oo_load_mh.oo_.mom,'convergence')
                    oo_.mom.convergence = oo_load_mh.oo_.mom.convergence;
                end
            end
        end
        % statistics and plots for posterior draws
        if options_mom_.mh_replic || (options_mom_.load_mh_file && ~options_mom_.load_results_after_load_mh)
            [~,oo_.mom] = marginal_density(M_, options_mom_, estim_params_, oo_.mom, bayestopt_, 'method_of_moments');
            oo_.mom = GetPosteriorParametersStatistics(estim_params_, M_, options_mom_, bayestopt_, oo_.mom, prior_dist_names);
            if ~options_mom_.nograph
                oo_.mom = PlotPosteriorDistributions(estim_params_, M_, options_mom_, bayestopt_, oo_.mom);
            end
            [oo_.mom.posterior.metropolis.mean,oo_.mom.posterior.metropolis.Variance] = GetPosteriorMeanVariance(options_mom_,M_);
        elseif options_mom_.load_mh_file && options_mom_.load_results_after_load_mh
            % load fields from previous MCMC run stored in results-file
            field_names={'posterior_mode','posterior_std_at_mode',...% fields set by marginal_density
                         'posterior_mean','posterior_hpdinf','posterior_hpdsup','posterior_median','posterior_variance','posterior_std','posterior_deciles','posterior_density',...% fields set by GetPosteriorParametersStatistics
                         'prior_density',...% fields set by PlotPosteriorDistributions
                        };
            for field_iter=1:size(field_names,2)
                if isfield(oo_load_mh.oo_.mom,field_names{1,field_iter})
                    oo_.mom.(field_names{1,field_iter}) = oo_load_mh.oo_.mom.(field_names{1,field_iter});
                end
            end
            if isfield(oo_load_mh.oo_.mom,'MarginalDensity') && isfield(oo_load_mh.oo_.mom.MarginalDensity,'ModifiedHarmonicMean') % field set by marginal_density            
                oo_.mom.MarginalDensity.ModifiedHarmonicMean = oo_load_mh.oo_.mom.MarginalDensity.ModifiedHarmonicMean;
            end            
            if isfield(oo_load_mh.oo_.mom,'posterior') && isfield(oo_load_mh.oo_.mom.posterior,'metropolis') % field set by GetPosteriorMeanVariance
                oo_.mom.posterior.metropolis = oo_load_mh.oo_.mom.posterior.metropolis;
            end
        end
        [error_flag,~,options_mom_]= metropolis_draw(1,options_mom_,estim_params_,M_);
        if ~(~isempty(options_mom_.sub_draws) && options_mom_.sub_draws==0)
            % THIS IS PROBABLY NOT USEFUL HERE AND CAN BE REMOVED (PREPROCESSOR: REMOVE bayesian_irf, moments_varendo)
            %if options_mom_.bayesian_irf
            %    if error_flag
            %        error('method_of_moments: Cannot compute the posterior IRFs!');
            %    end
            %    PosteriorIRF('posterior','method_of_moments::mcmc');
            %end
            % if options_mom_.moments_varendo
            %     if error_flag
            %         error('method_of_moments: Cannot compute the posterior moments for the endogenous variables!');
            %     end
            %     if options_mom_.load_mh_file && options_mom_.mh_replic==0 %user wants to recompute results
            %        [MetropolisFolder, info] = CheckPath('metropolis',options_mom_.dirname);
            %        if ~info
            %            generic_post_data_file_name={'Posterior2ndOrderMoments','decomposition','PosteriorVarianceDecomposition','correlation','PosteriorCorrelations','conditional decomposition','PosteriorConditionalVarianceDecomposition'};
            %            for ii=1:length(generic_post_data_file_name)
            %                delete_stale_file([MetropolisFolder filesep M_.fname '_' generic_post_data_file_name{1,ii} '*']);
            %            end
            %            % restore compatibility for loading pre-4.6.2 runs where estim_params_ was not saved; see 6e06acc7 and !1944
            %            NumberOfDrawsFiles = length(dir([M_.dname '/metropolis/' M_.fname '_posterior_draws*' ]));
            %            if NumberOfDrawsFiles>0
            %                temp=load([M_.dname '/metropolis/' M_.fname '_posterior_draws1']);
            %                if ~isfield(temp,'estim_params_')
            %                    for file_iter=1:NumberOfDrawsFiles
            %                        save([M_.dname '/metropolis/' M_.fname '_posterior_draws' num2str(file_iter)],'estim_params_','-append')
            %                    end
            %                end
            %            end
            %        end
            %     end
            %     oo_ = compute_moments_varendo('posterior',options_,M_,oo_,var_list_);
            % end            
        else
            fprintf('''sub_draws'' was set to 0. Skipping posterior computations.');
        end
        xparam1 = get_posterior_parameters('mean',M_,estim_params_,oo_.mom,options_);
    end
    % MAYBE USEFUL????
    % % Posterior correlations
    % extreme_corr_bound = 0.7;
    % if  ~isnan(extreme_corr_bound)
    %     tril_para_correlation_matrix=tril(para_correlation_matrix,-1);
    %     [row_index,col_index]=find(abs(tril_para_correlation_matrix)>extreme_corr_bound);
    %     extreme_corr_params=cell(length(row_index),3);
    %     for i=1:length(row_index)
    %         extreme_corr_params{i,1}=char(parameter_names(row_index(i),:));
    %         extreme_corr_params{i,2}=char(parameter_names(col_index(i),:));
    %         extreme_corr_params{i,3}=tril_para_correlation_matrix(row_index(i),col_index(i));
    %     end
    % end
    % disp(' ');
    % disp(['Correlations of Parameters (at Posterior Mode) > ',num2str(extreme_corr_bound)]);
    % disp(extreme_corr_params)
end


% -------------------------------------------------------------------------
% display final estimation results
% -------------------------------------------------------------------------
M_ = set_all_parameters(xparam1,estim_params_,M_); % update parameters
[~, ~, ~, ~, ~, oo_.mom.Q, oo_.mom.model_moments, oo_.mom.model_moments_params_derivs, oo_.mom.irf_model_varobs] = objective_function(xparam1, oo_.mom.data_moments, oo_.mom.weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);  % store final results in oo_.mom
if strcmp(options_mom_.mom.mom_method,'SMM') || strcmp(options_mom_.mom.mom_method,'GMM')
    % J test
    oo_.mom.J_test = mom.Jtest(xparam1, objective_function, oo_.mom.Q, oo_.mom.model_moments, oo_.mom.m_data, oo_.mom.data_moments, oo_.mom.weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
elseif strcmp(options_mom_.mom.mom_method,'IRF_MATCHING')
    if ~options_mom_.nograph
        mom.graph_comparison_irfs(M_.matched_irfs,oo_.mom.irf_model_varobs,options_mom_.varobs_id,options_mom_.irf,options_mom_.relative_irf,M_.endo_names,M_.endo_names_tex,M_.exo_names,M_.exo_names_tex,M_.dname,M_.fname,options_mom_.graph_format,options_mom_.TeX,options_mom_.nodisplay,options_mom_.figures.textwidth)
    end
end
% display comparison of model moments/IRFs and data moments/IRFs
mom.display_comparison_moments_irfs(M_, options_mom_, oo_.mom.data_moments, oo_.mom.model_moments);
% save results to _mom_results.mat
save([M_.dname filesep 'method_of_moments' filesep M_.fname '_mom_results.mat'], 'oo_', 'options_mom_', 'M_', 'estim_params_', 'bayestopt_');

fprintf('\n==== Method of Moments Estimation (%s) Completed ====\n\n',options_mom_.mom.mom_method);

% -------------------------------------------------------------------------
% clean up
% -------------------------------------------------------------------------
warning(orig_warning_state); %reset warning state
if isoctave && isfield(options_mom_, 'prior_restrictions') && ...
   isfield(options_mom_.prior_restrictions, 'routine')
    % Octave crashes if it tries to save function handles (to the _results.mat file)
    % See https://savannah.gnu.org/bugs/?43215
    options_mom_.prior_restrictions.routine = [];
end
if strcmp(options_mom_.mom.mom_method,'SMM') || strcmp(options_mom_.mom.mom_method,'GMM')
    if isfield(oo_.mom,'irf_model_varobs') && isempty(oo_.mom.irf_model_varobs)
        oo_.mom = rmfield(oo_.mom,'irf_model_varobs'); % remove empty field
    end
end
if strcmp(options_mom_.mom.mom_method,'IRF_MATCHING') && ~isempty(options_mom_.mom.irf_matching_file.path) && ~strcmp(options_mom_.mom.irf_matching_file.path,'.')
    rmpath(options_mom_.irf_matching_file.path); % remove path to irf_matching_file
end