function [oo_, options_mom_, M_] = run(bayestopt_, options_, oo_, estim_params_, M_, options_mom_)
% function [oo_, options_mom_, M_] = run(bayestopt_, options_, oo_, estim_params_, M_, options_mom_)
% -------------------------------------------------------------------------
% This function performs a method of moments estimation with the following steps:
%  o Checking if required structures and options exist
%  o Preparing local options_mom_ structure
%  o Checking the options and the compatibility of the settings
%  o Initializations of variables, orderings and state space representation
%  o Checks and transformations for matched moments structure
%  o Checks and transformations for estimated parameters, priors, and bounds
%  o Checks and transformations for data
%  o Checks for objective function at initial parameters
%  o GMM/SMM: iterated method of moments estimation
%  o GMM/SMM: J-Test and fit of moments% 
%  o Display of results
%  o Clean up
% -------------------------------------------------------------------------
% This function is inspired by replication codes accompanied to the following papers:
% GMM/SMM:
%  o Andreasen, Fernández-Villaverde, Rubio-Ramírez (2018): "The Pruned State-Space System for Non-Linear DSGE Models: Theory and Empirical Applications", Review of Economic Studies, 85(1):1-49.
%  o Born, Pfeifer (2014): "Risk Matters: Comment", American Economic Review, 104(12):4231-4239.
%  o Mutschler (2018): "Higher-order statistics for DSGE models", Econometrics and Statistics, 6:44-56.
% =========================================================================
% INPUTS
%  o bayestopt_:     [structure] information about priors
%  o options_:       [structure] information about global options
%  o oo_:            [structure] storage for results
%  o estim_params_:  [structure] information about estimated parameters
%  o M_:             [structure] information about model with
%                     o matched_moments:  [cell] information about selected moments to match in GMM/SMM estimation
%                                                vars: matched_moments{:,1});
%                                                lead/lags: matched_moments{:,2};
%                                                powers: matched_moments{:,3};
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
%  o cellofchararraymaxlength
%  o check_for_calibrated_covariances
%  o check_prior_bounds
%  o check_prior_stderr_corr
%  o check_steady_state_changes_parameters
%  o check_varobs_are_endo_and_declared_once
%  o display_estimation_results_table
%  o do_parameter_initialization
%  o dyn_latex_table
%  o dynare_minimize_objective
%  o dyntable
%  o get_all_parameters
%  o get_dynare_random_generator_state
%  o get_matrix_entries_for_psd_check
%  o M_.fname '_prior_restrictions'
%  o makedataset
%  o mom.check_plot
%  o mom.default_option_mom_values
%  o mom.get_data_moments
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
%  o warning_config
% =========================================================================
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
% -------------------------------------------------------------------------
% Maintaining Author(s):
% o Willi Mutschler (willi@mutschler.eu)
% o Johannes Pfeifer (johannes.pfeifer@unibw.de)
% =========================================================================

% -------------------------------------------------------------------------
% TO DO LISTS
% -------------------------------------------------------------------------
% GENERAL
% - document all options in manual
% - document analytic_jacobian better
% - make endogenous_prior_restrictions work
% - dirname option to save output to different directory not yet implemented
% - create test for prior restrictions file
% - add mode_file option
% - implement penalty objective
% - test optimizers
% GMM/SMM
% - speed up pruned_state_space_system (by using doubling with old initial values, hardcoding zeros, other "tricks" used in e.g. nlma)
% - add option to use autocorrelations (we have useautocorr in identification toolbox already)
% - SMM with extended path
% - deal with measurement errors (once @wmutschl has implemented this in identification toolbox)
% - display scaled moments
% - enable first moments despite prefilter
% - do "true" Bayesian GMM and SMM not only penalized

fprintf('\n==== Method of Moments Estimation (%s) ====\n\n',options_mom_.mom.mom_method)


% -------------------------------------------------------------------------
% checks if required structures exist
% -------------------------------------------------------------------------
if isempty(estim_params_) % structure storing the info about estimated parameters in the estimated_params block
    if ~(isfield(estim_params_,'nvx') && (size(estim_params_.var_exo,1)+size(estim_params_.var_endo,1)+size(estim_params_.corrx,1)+size(estim_params_.corrn,1)+size(estim_params_.param_vals,1))==0)
        error('method_of_moments: You need to provide an ''estimated_params'' block!')
    else
        error('method_of_moments: The ''estimated_params'' block must not be empty!')
    end
end
if strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')
    if ~isfield(M_,'matched_moments') || isempty(M_.matched_moments) % structure storing the moments used for GMM and SMM estimation
        error('method_of_moments: You need to provide a ''matched_moments'' block for ''mom_method=%s''!',options_mom_.mom.mom_method)
    end
end
if (~isempty(estim_params_.var_endo) || ~isempty(estim_params_.corrn)) && strcmp(options_mom_.mom.mom_method, 'GMM')
    error('method_of_moments: GMM estimation does not support measurement error(s) yet. Please specifiy them as a structural shock!')
end
doBayesianEstimation = [estim_params_.var_exo(:,5); estim_params_.var_endo(:,5); estim_params_.corrx(:,6); estim_params_.corrn(:,6); estim_params_.param_vals(:,5)];
if all(doBayesianEstimation~=0)
    doBayesianEstimation = true;
elseif all(doBayesianEstimation==0)
    doBayesianEstimation = false;
else
    error('method_of_moments: Estimation must be either fully Frequentist or fully Bayesian. Maybe you forgot to specify a prior distribution!')
end
if ~isfield(options_,'varobs')
    error('method_of_moments: VAROBS statement is missing!')
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
options_mom_ = mom.default_option_mom_values(options_mom_, options_, M_.dname, doBayesianEstimation);


% -------------------------------------------------------------------------
% workarounds
% -------------------------------------------------------------------------
% The TeX option crashes MATLAB R2014a run with "-nodisplay" option
% (as is done from the testsuite).
% Since we can’t directly test whether "-nodisplay" has been passed,
% we test for the "TOP_TEST_DIR" environment variable, which is set
% by the testsuite.
% Note that it was not tested whether the crash happens with more
% recent MATLAB versions, so when OLD_MATLAB_VERSION is increased,
% one should make a test before removing this workaround.
if options_.TeX && ~isoctave && matlab_ver_less_than('8.4') && ~isempty(getenv('TOP_TEST_DIR'))
    warning('Disabling TeX option due to a bug in MATLAB R2014a with -nodisplay')
    options_.TeX = false;
end
if isfield(options_mom_, 'TeX') && options_mom_.TeX && ~isoctave && matlab_ver_less_than('8.4') && ~isempty(getenv('TOP_TEST_DIR'))
    warning('Disabling TeX option due to a bug in MATLAB R2014a with -nodisplay')
    options_mom_.TeX = false;
end
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
% checks on settings
% -------------------------------------------------------------------------
if strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')
    if numel(options_mom_.nobs) > 1
        error('method_of_moments: Recursive estimation and forecast for samples is not supported. Please set an integer as ''nobs''!');
    end
    if numel(options_mom_.first_obs) > 1
        error('method_of_moments: Recursive estimation and forecast for samples is not supported. Please set an integer as ''first_obs''!');
    end
end
if options_mom_.order < 1
    error('method_of_moments: The order of the Taylor approximation cannot be 0!')
end
if options_mom_.order > 2
    fprintf('Dynare will use ''k_order_solver'' as the order>2\n');
    options_mom_.k_order_solver = true;
end
if strcmp(options_mom_.mom.mom_method,'SMM')
    if options_mom_.mom.simulation_multiple < 1
        fprintf('The simulation horizon is shorter than the data. Dynare resets the simulation_multiple to 7.\n')
        options_mom_.mom.simulation_multiple = 7;
    end
end
if strcmp(options_mom_.mom.mom_method,'GMM')
    % require pruning with GMM at higher order
    if options_mom_.order > 1 && ~options_mom_.pruning
        fprintf('GMM at higher order only works with pruning, so we set pruning option to 1.\n');
        options_mom_.pruning = true;
    end
    if options_mom_.order > 3
        error('method_of_moments: Perturbation orders higher than 3 are not implemented for GMM estimation, try using SMM!');
    end
end
if options_mom_.mom.analytic_jacobian && ~strcmp(options_mom_.mom.mom_method,'GMM')
    options_mom_.mom.analytic_jacobian = false;
    fprintf('\n''analytic_jacobian'' option will be dismissed as it only works with GMM.\n');
end


% -------------------------------------------------------------------------
% initializations
% -------------------------------------------------------------------------
% create output directories to store results
CheckPath('method_of_moments',M_.dname);
CheckPath('graphs',options_mom_.dirname);
% initialize options that might change
options_mom_.mom.compute_derivs = false; % flag to compute derivs in objective function (might change for GMM with either analytic_standard_errors or analytic_jacobian (dependent on optimizer))
options_mom_.mom.vector_output = false;  % specifies whether the objective function returns a vector
% decision rule
oo_.dr = set_state_space(oo_.dr,M_,options_mom_); % get state-space representation
oo_.mom.obs_var = []; % create index of observed variables in DR order
for i = 1:options_mom_.obs_nbr
    oo_.mom.obs_var = [oo_.mom.obs_var; find(strcmp(options_mom_.varobs{i}, M_.endo_names(oo_.dr.order_var)))];
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
    not_observed_variables=setdiff(oo_.dr.inv_order_var([M_.matched_moments{:,1}]),oo_.mom.obs_var);
    if ~isempty(not_observed_variables)
        skipline;
        error('method_of_moments: You specified moments involving %s, but it is not a varobs!',M_.endo_names{oo_.dr.order_var(not_observed_variables)})
    end
end


% -------------------------------------------------------------------------
% estimated parameters: checks and transformations on values, priors, bounds
% -------------------------------------------------------------------------
% set priors and bounds over the estimated parameters
[xparam0, estim_params_, bayestopt_, lb, ub, M_] = set_prior(estim_params_, M_, options_mom_);
number_of_estimated_parameters = length(xparam0);
hessian_xparam0 = []; % initialize hessian

% check if enough moments for estimation
if strcmp(options_mom_.mom.mom_method, 'GMM') || strcmp(options_mom_.mom.mom_method, 'SMM')
    if options_mom_.mom.mom_nbr < length(xparam0)
        skipline;
        error('method_of_moments: There must be at least as many moments as parameters for a %s estimation!',options_mom_.mom.mom_method);
    end
    skipline(2);
end

% check if a _prior_restrictions.m file exists
if exist([M_.fname '_prior_restrictions.m'],'file')
    options_mom_.prior_restrictions.status = 1;
    options_mom_.prior_restrictions.routine = str2func([M_.fname '_prior_restrictions']);
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
estim_params_ = check_for_calibrated_covariances(xparam0,estim_params_,M_);

% checks on parameter calibration and initialization
xparam_calib = get_all_parameters(estim_params_,M_); % get calibrated parameters
if ~any(isnan(xparam_calib)) % all estimated parameters are calibrated
    estim_params_.full_calibration_detected = true;
else
    estim_params_.full_calibration_detected = false;
end
if options_mom_.use_calibration_initialization % set calibration as starting values
    if ~isempty(bayestopt_) && ~doBayesianEstimation && any(all(isnan([xparam_calib xparam0]),2))
        error('method_of_moments: When using the use_calibration option with %s without prior, the parameters must be explicitly initialized!',options_mom_.mom.mom_method);
    else
        [xparam0,estim_params_] = do_parameter_initialization(estim_params_,xparam_calib,xparam0); % get explicitly initialized parameters that have precedence over calibrated values
    end
end

% check initialization
if ~isempty(bayestopt_) && ~doBayesianEstimation && any(isnan(xparam0))
    error('method_of_moments: Frequentist %s requires all estimated parameters to be initialized, either in an estimated_params or estimated_params_init-block!',options_mom_.mom.mom_method);
end

% set and check parameter bounds
if ~isempty(bayestopt_) && doBayesianEstimation
    % plot prior densities
    if ~options_mom_.nograph && options_mom_.plot_priors
        if strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')
            plot_priors(bayestopt_orig,M_,estim_params_,options_mom_,'Original priors'); % only for visual inspection (not saved to disk, because overwritten in next call to plot_priors)
            plot_priors(bayestopt_,M_,estim_params_,options_mom_,'Laplace approximated priors');
            clear('bayestopt_orig'); % make sure stale structure cannot be used
        end
    end
    % set prior bounds
    Bounds = prior_bounds(bayestopt_, options_mom_.prior_trunc);
    Bounds.lb = max(Bounds.lb,lb);
    Bounds.ub = min(Bounds.ub,ub);
else
    % no priors are declared so Dynare will estimate the parameters with
    % Frequentist methods using inequality constraints for the parameters
    Bounds.lb = lb;
    Bounds.ub = ub;
    if options_mom_.mom.penalized_estimator
        fprintf('Penalized estimation turned off as you did not declare priors\n')
        options_mom_.mom.penalized_estimator = 0;
    end
end

% set correct bounds for standard deviations and correlations
Bounds = mom.set_correct_bounds_for_stderr_corr(estim_params_,Bounds);

% test if initial values of the estimated parameters are all between the prior lower and upper bounds
if options_mom_.use_calibration_initialization
    try
        check_prior_bounds(xparam0,Bounds,M_,estim_params_,options_mom_,bayestopt_);
    catch last_error
        fprintf('Cannot use parameter values from calibration as they violate the prior bounds.')
        rethrow(last_error);
    end
else
    check_prior_bounds(xparam0,Bounds,M_,estim_params_,options_mom_,bayestopt_);
end

% check for positive definiteness
estim_params_ = get_matrix_entries_for_psd_check(M_,estim_params_);

% set sigma_e_is_diagonal flag (needed if the shocks block is not declared in the mod file)
M_.sigma_e_is_diagonal = true;
if estim_params_.ncx || any(nnz(tril(M_.Correlation_matrix,-1))) || isfield(estim_params_,'calibrated_covariances')
    M_.sigma_e_is_diagonal = false;
end

% storing prior parameters in results
oo_.mom.prior.mean = bayestopt_.p1;
oo_.mom.prior.mode = bayestopt_.p5;
oo_.mom.prior.variance = diag(bayestopt_.p2.^2);
oo_.mom.prior.hyperparameters.first = bayestopt_.p6;
oo_.mom.prior.hyperparameters.second = bayestopt_.p7;

% set all parameters
M_ = set_all_parameters(xparam0,estim_params_,M_);

% provide warning if there is NaN in parameters
test_for_deep_parameters_calibration(M_);

if doBayesianEstimation
    % warning if prior allows that stderr parameters are negative or corr parameters are outside the unit circle
    check_prior_stderr_corr(estim_params_,bayestopt_);

    % check value of prior density
    [~,~,~,info]= priordens(xparam0,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4);
    if any(info)
        fprintf('The prior density evaluated at the initial values is Inf for the following parameters: %s\n',bayestopt_.name{info,1})
        error('The initial value of the prior is -Inf!')
    end
end


% -------------------------------------------------------------------------
% datafile: checks and transformations
% -------------------------------------------------------------------------
% Build dataset
if strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')
    % Check if datafile has same name as mod file
    [~,name,~] = fileparts(options_mom_.datafile);
    if strcmp(name,M_.fname)
        error('method_of_moments: ''datafile'' and mod file are not allowed to have the same name; change the name of the ''datafile''!')
    end
    dataset_ = makedataset(options_mom_);
    % set options for old interface from the ones for new interface
    if ~isempty(dataset_)
        options_mom_.nobs = dataset_.nobs;
    end
    % Check length of data for estimation of second moments
    if options_mom_.ar > options_mom_.nobs+1
        error('method_of_moments: Dataset is too short to compute higher than first moments!');
    end
    % Provide info on data moments handling
    fprintf('Computing data moments. Note that NaN values in the moments (due to leads and lags or missing data) are replaced by the mean of the corresponding moment.\n');
    % Get data moments for the method of moments
    [oo_.mom.data_moments, oo_.mom.m_data] = mom.get_data_moments(dataset_.data, oo_, M_.matched_moments, options_mom_);
    if ~isreal(dataset_.data)
        error('method_of_moments: The data moments contain complex values!')
    end
end


% -------------------------------------------------------------------------
% SMM: Get shock series fand set variance correction factor
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
% Step 5: checks for steady state at initial parameters
% -------------------------------------------------------------------------

% setting steadystate_check_flag option
if options_mom_.steadystate.nocheck
    steadystate_check_flag = 0;
else
    steadystate_check_flag = 1;
end

old_steady_params=M_.params; %save initial parameters for check if steady state changes param values
% Check steady state at initial model parameter values
[oo_.steady_state, new_steady_params, info] = evaluate_steady_state(oo_.steady_state,[oo_.exo_steady_state; oo_.exo_det_steady_state],M_,options_mom_,steadystate_check_flag);
if info(1)
    fprintf('\nmethod_of_moments: The steady state at the initial parameters cannot be computed.\n')
    print_info(info, 0, options_mom_);
end

% check whether steady state file changes estimated parameters
if isfield(estim_params_,'param_vals') && ~isempty(estim_params_.param_vals)
    Model_par_varied=M_; %store M_ structure
    
    Model_par_varied.params(estim_params_.param_vals(:,1))=Model_par_varied.params(estim_params_.param_vals(:,1))*1.01; %vary parameters
    [~, new_steady_params_2] = evaluate_steady_state(oo_.steady_state,[oo_.exo_steady_state; oo_.exo_det_steady_state],Model_par_varied,options_mom_,true);
    
    changed_par_indices=find((old_steady_params(estim_params_.param_vals(:,1))-new_steady_params(estim_params_.param_vals(:,1))) ...
        | (Model_par_varied.params(estim_params_.param_vals(:,1))-new_steady_params_2(estim_params_.param_vals(:,1))));
    
    if ~isempty(changed_par_indices)
        fprintf('\nThe steady state file internally changed the values of the following estimated parameters:\n')
        disp(char(M_.param_names(estim_params_.param_vals(changed_par_indices,1))))
        fprintf('This will override parameter values and may lead to wrong results.\n')
        fprintf('Check whether this is really intended.\n')
        warning('The steady state file internally changes the values of the estimated parameters.')
        if strcmp(options_mom_.mom.mom_method,'GMM') && options_mom_.mom.analytic_standard_errors
            fprintf('For analytical standard errors, the parameter-Jacobians of the dynamic model and of the steady-state will be computed numerically\n'),
            fprintf('(re-set options_mom_.analytic_derivation_mode= -2)'),
            options_mom_.analytic_derivation_mode= -2;
        end
    end
end

% display warning if some parameters are still NaN
test_for_deep_parameters_calibration(M_);

% If steady state of observed variables is non zero, set noconstant equal 0
if all(abs(oo_.steady_state(oo_.dr.order_var(oo_.mom.obs_var)))<1e-9)
    options_mom_.noconstant = 0; %identifying the constant based on just the initial parameter value is not feasible
else
    options_mom_.noconstant = 0;
end

% -------------------------------------------------------------------------
% Step 6: checks for objective function at initial parameters
% -------------------------------------------------------------------------
objective_function = str2func('mom.objective_function');

try
    % Check for NaN or complex values of moment-distance-funtion evaluated
    % at initial parameters and identity weighting matrix    
    oo_.mom.Sw = eye(options_mom_.mom.mom_nbr);
    tic_id = tic;    
    [fval, info, ~, ~, ~, oo_, M_] = feval(objective_function, xparam0, Bounds, oo_, estim_params_, M_, options_mom_);
    elapsed_time = toc(tic_id);    
    if isnan(fval)
        error('method_of_moments: The initial value of the objective function is NaN')
    elseif imag(fval)
        error('method_of_moments: The initial value of the objective function is complex')
    end
    if info(1) > 0
        disp('method_of_moments: Error in computing the objective function for initial parameter values')
        print_info(info, options_mom_.noprint, options_mom_)
    end
    fprintf('Initial value of the moment objective function with %4.1f times identity weighting matrix: %6.4f \n\n', options_mom_.mom.weighting_matrix_scaling_factor, fval);
    fprintf('Time required to compute objective function once: %5.4f seconds \n', elapsed_time);
    
catch last_error% if check fails, provide info on using calibration if present
    if estim_params_.full_calibration_detected %calibrated model present and no explicit starting values
        skipline(1);
        fprintf('There was an error in computing the moments for initial parameter values.\n')
        fprintf('If this is not a problem with the setting of options (check the error message below),\n')
        fprintf('you should try using the calibrated version of the model as starting values. To do\n')
        fprintf('this, add an empty estimated_params_init-block with use_calibration option immediately before the estimation\n')
        fprintf('command (and after the estimated_params-block so that it does not get overwritten):\n');
        skipline(2);
    end
    rethrow(last_error);
end

% -------------------------------------------------------------------------
% Step 7a: Method of moments estimation: print some info
% -------------------------------------------------------------------------
fprintf('\n---------------------------------------------------\n')
if strcmp(options_mom_.mom.mom_method,'SMM')
    fprintf('Simulated method of moments with');
elseif strcmp(options_mom_.mom.mom_method,'GMM')
    fprintf('General method of moments with');
end
if options_mom_.prefilter
    fprintf('\n  - centered moments (prefilter=1)');
else
    fprintf('\n  - uncentered moments (prefilter=0)');
end
if options_mom_.mom.penalized_estimator
    fprintf('\n  - penalized estimation using deviation from prior mean and weighted with prior precision');
end

for i = 1:length(options_mom_.optimizer_vec)
    if i == 1
        str = '- optimizer (mode_compute';
    else
        str = '            (additional_optimizer_steps';
    end
    switch options_mom_.optimizer_vec{i}
        case 0
            fprintf('\n  %s=0): no minimization',str);
        case 1
            fprintf('\n  %s=1): fmincon',str);
        case 2
            fprintf('\n  %s=2): continuous simulated annealing',str);
        case 3
            fprintf('\n  %s=3): fminunc',str);
        case 4
            fprintf('\n  %s=4): csminwel',str);
        case 5
            fprintf('\n  %s=5): newrat',str);
        case 6
            fprintf('\n  %s=6): gmhmaxlik',str);
        case 7
            fprintf('\n  %s=7): fminsearch',str);
        case 8
            fprintf('\n  %s=8): Dynare Nelder-Mead simplex',str);
        case 9
            fprintf('\n  %s=9): CMA-ES',str);
        case 10
            fprintf('\n  %s=10): simpsa',str);
        case 11
            error('\nmethod_of_moments: online_auxiliary_filter (mode_compute=11) is only supported with likelihood-based estimation techniques');
        case 12
            fprintf('\n  %s=12): particleswarm',str);
        case 101
            fprintf('\n  %s=101): SolveOpt',str);
        case 102
            fprintf('\n  %s=102): simulannealbnd',str);
        case 13
            fprintf('\n  %s=13): lsqnonlin',str);
        otherwise
            if ischar(options_mom_.optimizer_vec{i})
                fprintf('\n  %s=%s): user-defined',str,options_mom_.optimizer_vec{i});
            else
                error('method_of_moments: Unknown optimizer, please contact the developers ')
            end
    end
    if options_mom_.silent_optimizer
        fprintf(' (silent)');
    end
    if strcmp(options_mom_.mom.mom_method,'GMM') && options_mom_.mom.analytic_jacobian && ismember(options_mom_.optimizer_vec{i},options_mom_.mom.analytic_jacobian_optimizers)
        fprintf(' (using analytical Jacobian)');
    end
end
fprintf('\n  - perturbation order:        %d', options_mom_.order)
if options_mom_.order > 1 && options_mom_.pruning
    fprintf(' (with pruning)')
end
if strcmp(options_mom_.mom.mom_method,'GMM') && options_mom_.mom.analytic_standard_errors
    fprintf('\n  - standard errors:           analytic derivatives');
else
    fprintf('\n  - standard errors:           numerical derivatives');
end
fprintf('\n  - number of matched moments: %d', options_mom_.mom.mom_nbr);
fprintf('\n  - number of parameters:      %d', length(xparam0));
fprintf('\n\n');

% -------------------------------------------------------------------------
% Step 7b: Iterated method of moments estimation
% -------------------------------------------------------------------------
if size(options_mom_.mom.weighting_matrix,1)>1 && ~(any(strcmpi('diagonal',options_mom_.mom.weighting_matrix)) || any(strcmpi('optimal',options_mom_.mom.weighting_matrix)))
    fprintf('\nYou did not specify the use of an optimal or diagonal weighting matrix. There is no point in running an iterated method of moments.\n')
end

for stage_iter=1:size(options_mom_.mom.weighting_matrix,1)
    fprintf('Estimation stage %u\n',stage_iter);
    Woptflag = false;
    switch lower(options_mom_.mom.weighting_matrix{stage_iter})
        case 'identity_matrix'
            fprintf('  - identity weighting matrix\n');
            weighting_matrix = eye(options_mom_.mom.mom_nbr);            
        case 'diagonal'
            fprintf('  - diagonal of optimal weighting matrix (Bartlett kernel with %d lags)\n', options_mom_.mom.bartlett_kernel_lag);
            if stage_iter == 1
                fprintf('    and using data-moments as initial estimate of model-moments\n');
                weighting_matrix = diag(diag(  mom.optimal_weighting_matrix(oo_.mom.m_data, oo_.mom.data_moments, options_mom_.mom.bartlett_kernel_lag)  ));
            else
                fprintf('    and using previous stage estimate of model-moments\n');
                weighting_matrix = diag(diag(  mom.optimal_weighting_matrix(oo_.mom.m_data, oo_.mom.model_moments, options_mom_.mom.bartlett_kernel_lag)  ));
            end            
        case 'optimal'
            fprintf('  - optimal weighting matrix (Bartlett kernel with %d lags)\n', options_mom_.mom.bartlett_kernel_lag);
            if stage_iter == 1
                fprintf('    and using data-moments as initial estimate of model-moments\n');
                weighting_matrix = mom.optimal_weighting_matrix(oo_.mom.m_data, oo_.mom.data_moments, options_mom_.mom.bartlett_kernel_lag);
            else
                fprintf('    and using previous stage estimate of model-moments\n');
                weighting_matrix = mom.optimal_weighting_matrix(oo_.mom.m_data, oo_.mom.model_moments, options_mom_.mom.bartlett_kernel_lag);
                Woptflag = true;
            end            
        otherwise %user specified matrix in file
            fprintf('  - user-specified weighting matrix\n');
            try
                load(options_mom_.mom.weighting_matrix{stage_iter},'weighting_matrix')
            catch
                error(['method_of_moments: No matrix named ''weighting_matrix'' could be found in ',options_mom_.mom.weighting_matrix{stage_iter},'.mat'])
            end
            [nrow, ncol] = size(weighting_matrix);
            if ~isequal(nrow,ncol) || ~isequal(nrow,length(oo_.mom.data_moments)) %check if square and right size
                error(['method_of_moments: weighting_matrix must be square and have ',num2str(length(oo_.mom.data_moments)),' rows and columns'])
            end            
    end
    try %check for positive definiteness of weighting_matrix
        oo_.mom.Sw = chol(weighting_matrix);
    catch
        error('method_of_moments: Specified weighting_matrix is not positive definite. Check whether your model implies stochastic singularity.')
    end

    for optim_iter= 1:length(options_mom_.optimizer_vec)
        options_mom_.current_optimizer = options_mom_.optimizer_vec{optim_iter};
        if options_mom_.optimizer_vec{optim_iter}==0
            xparam1=xparam0; %no minimization, evaluate objective at current values
            fval = feval(objective_function, xparam1, Bounds, oo_, estim_params_, M_, options_mom_);
        else
            if options_mom_.optimizer_vec{optim_iter}==13
                options_mom_.mom.vector_output = true;                
            else
                options_mom_.mom.vector_output = false;                
            end
            if strcmp(options_mom_.mom.mom_method,'GMM') && options_mom_.mom.analytic_jacobian && ismember(options_mom_.optimizer_vec{optim_iter},options_mom_.mom.analytic_jacobian_optimizers) %do this only for gradient-based optimizers
                options_mom_.mom.compute_derivs = true;
            else
                options_mom_.mom.compute_derivs = false;
            end
            
            [xparam1, fval, exitflag] = dynare_minimize_objective(objective_function, xparam0, options_mom_.optimizer_vec{optim_iter}, options_mom_, [Bounds.lb Bounds.ub], bayestopt_.name, bayestopt_, [],...
                                                                  Bounds, oo_, estim_params_, M_, options_mom_);
            if options_mom_.mom.vector_output
                fval = fval'*fval;
            end
        end
        fprintf('\nStage %d Iteration %d: value of minimized moment distance objective function: %12.10f.\n',stage_iter,optim_iter,fval)
        if options_mom_.mom.verbose
            oo_.mom=display_estimation_results_table(xparam1,NaN(size(xparam1)),M_,options_mom_,estim_params_,bayestopt_,oo_.mom,prior_dist_names,sprintf('%s (STAGE %d ITERATION %d) VERBOSE',options_mom_.mom.mom_method,stage_iter,optim_iter),sprintf('verbose_%s_stage_%d_iter_%d',lower(options_mom_.mom.mom_method),stage_iter,optim_iter));
        end
        xparam0=xparam1;
    end
    options_mom_.mom.vector_output = false;    
    % Update M_ and DynareResults (in particular to get oo_.mom.model_moments)    
    M_ = set_all_parameters(xparam1,estim_params_,M_);
    if strcmp(options_mom_.mom.mom_method,'GMM') && options_mom_.mom.analytic_standard_errors
        options_mom_.mom.compute_derivs = true; % for GMM we compute derivatives analytically in the objective function with this flag        
    end
    [fval, ~, ~,~,~, oo_] = feval(objective_function, xparam1, Bounds, oo_, estim_params_, M_, options_mom_);
    options_mom_.mom.compute_derivs = false; % reset to not compute derivatives in objective function during optimization
    
    SE = mom.standard_errors(xparam1, objective_function, Bounds, oo_, estim_params_, M_, options_mom_, Woptflag);
    
    % Store results in output structure
    oo_.mom = display_estimation_results_table(xparam1,SE,M_,options_mom_,estim_params_,bayestopt_,oo_.mom,prior_dist_names,sprintf('%s (STAGE %u)',options_mom_.mom.mom_method,stage_iter),sprintf('%s_stage_%u',lower(options_mom_.mom.mom_method),stage_iter));
end

% -------------------------------------------------------------------------
% Step 8: J test
% -------------------------------------------------------------------------
if options_mom_.mom.mom_nbr > length(xparam1)
    %get optimal weighting matrix for J test, if necessary
    if ~Woptflag
        W_opt = mom.optimal_weighting_matrix(oo_.mom.m_data, oo_.mom.model_moments, options_mom_.mom.bartlett_kernel_lag);
        oo_j=oo_;
        oo_j.mom.Sw = chol(W_opt);
        [fval] = feval(objective_function, xparam1, Bounds, oo_j, estim_params_, M_, options_mom_);
    end

    % Compute J statistic
    if strcmp(options_mom_.mom.mom_method,'SMM')    
        Variance_correction_factor = options_mom_.mom.variance_correction_factor;
    elseif strcmp(options_mom_.mom.mom_method,'GMM')
        Variance_correction_factor=1;
    end
    oo_.mom.J_test.j_stat          = dataset_.nobs*Variance_correction_factor*fval/options_mom_.mom.weighting_matrix_scaling_factor;
    oo_.mom.J_test.degrees_freedom = length(oo_.mom.model_moments)-length(xparam1);
    oo_.mom.J_test.p_val           = 1-chi2cdf(oo_.mom.J_test.j_stat, oo_.mom.J_test.degrees_freedom);
    fprintf('\nvalue of J-test statistic: %f\n',oo_.mom.J_test.j_stat)
    fprintf('p-value of J-test statistic: %f\n',oo_.mom.J_test.p_val)
end

% -------------------------------------------------------------------------
% Step 9: Display estimation results
% -------------------------------------------------------------------------
title = ['Comparison of data moments and model moments (',options_mom_.mom.mom_method,')'];
headers = {'Moment','Data','Model'};
for jm = 1:size(M_.matched_moments,1)
    lables_tmp = 'E[';
    lables_tmp_tex = 'E \left[ ';
    for jvar = 1:length(M_.matched_moments{jm,1})
        lables_tmp = [lables_tmp M_.endo_names{M_.matched_moments{jm,1}(jvar)}];
        lables_tmp_tex = [lables_tmp_tex, '{', M_.endo_names_tex{M_.matched_moments{jm,1}(jvar)}, '}'];
        if M_.matched_moments{jm,2}(jvar) ~= 0
            lables_tmp = [lables_tmp, '(', num2str(M_.matched_moments{jm,2}(jvar)), ')'];
            lables_tmp_tex = [lables_tmp_tex, '_{t', num2str(M_.matched_moments{jm,2}(jvar)), '}'];
        else
            lables_tmp_tex = [lables_tmp_tex, '_{t}'];
        end
        if M_.matched_moments{jm,3}(jvar) > 1
            lables_tmp = [lables_tmp, '^', num2str(M_.matched_moments{jm,3}(jvar))];
            lables_tmp_tex = [lables_tmp_tex, '^{', num2str(M_.matched_moments{jm,3}(jvar)) '}'];
        end
        if jvar == length(M_.matched_moments{jm,1})
            lables_tmp = [lables_tmp, ']'];
            lables_tmp_tex = [lables_tmp_tex, ' \right]'];
        else
            lables_tmp = [lables_tmp, '*'];
            lables_tmp_tex = [lables_tmp_tex, ' \times '];
        end
    end
    labels{jm,1} = lables_tmp;
    labels_TeX{jm,1} = lables_tmp_tex;
end
data_mat=[oo_.mom.data_moments oo_.mom.model_moments ];
dyntable(options_mom_, title, headers, labels, data_mat, cellofchararraymaxlength(labels)+2, 10, 7);
if options_mom_.TeX
    dyn_latex_table(M_, options_mom_, title, ['comparison_moments_', options_mom_.mom.mom_method], headers, labels_TeX, data_mat, cellofchararraymaxlength(labels)+2, 10, 7);
end

if options_mom_.mode_check.status
    mom.check_plot(objective_function,xparam1,SE,options_mom_,M_,estim_params_,Bounds,bayestopt_,...
        Bounds, oo_, estim_params_, M_, options_mom_)
end

fprintf('\n==== Method of Moments Estimation (%s) Completed ====\n\n',options_mom_.mom.mom_method)

% -------------------------------------------------------------------------
% Step 9: Clean up
% -------------------------------------------------------------------------
%reset warning state
warning_config;

if isoctave && isfield(options_, 'prior_restrictions') && ...
   isfield(options_.prior_restrictions, 'routine')
    % Octave crashes if it tries to save function handles (to the _results.mat file)
    % See https://savannah.gnu.org/bugs/?43215
    options_.prior_restrictions.routine = [];
end
