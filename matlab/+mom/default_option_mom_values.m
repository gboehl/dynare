function options_mom_ = default_option_mom_values(options_mom_, options_, dname, do_bayesian_estimation)
% options_mom_ = default_option_mom_values(options_mom_, options_, dname, do_bayesian_estimation)
% -------------------------------------------------------------------------
% Returns structure containing the options for method_of_moments command.
% Note 1: options_mom_ is local and contains default and user-specified
% values for all settings needed for the method of moments estimation.
% Some options, though, are set by the preprocessor into options_ and we
% copy these over. The idea is to be independent of options_ and have full
% control of the estimation instead of possibly having to deal with options
% chosen somewhere else in the mod file.
% Note 2: we call a "mode" the minimum of the objective function, i.e.
% the parameter vector that minimizes the distance between the moments/irfs
% computed from the model and the moments/irfs computed from the data.
% -------------------------------------------------------------------------
% INPUTS
%  o options_mom_:           [structure]  all user-specified settings (from the method_of_moments command)
%  o options_:               [structure]  global options
%  o dname:                  [string]     default name of directory to store results
%  o do_bayesian_estimation  [boolean]    indicator whether we do Bayesian estimation
% -------------------------------------------------------------------------
% OUTPUTS
%  o options_mom_:         [structure]  all user-specified and updated settings required for method_of_moments estimation
% -------------------------------------------------------------------------
% This function is called by
%  o mom.run
% -------------------------------------------------------------------------
% This function calls
%   o set_default_option
%   o user_has_matlab_license
%   o user_has_octave_forge_package
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


mom_method = options_mom_.mom.mom_method; % this is a required option


% -------------------------------------------------------------------------
% LIMITATIONS
% -------------------------------------------------------------------------

if options_.logged_steady_state || options_.loglinear
    error('method_of_moments: The loglinear option is not supported. Please append the required logged variables as auxiliary equations.');
else
    options_mom_.logged_steady_state = 0;
    options_mom_.loglinear = false;
end
if isfield(options_mom_,'hessian') && options_mom_.hessian.use_penalized_objective
    warning('method_of_moments: The ''use_penalized_objective_for_hessian'' option is not supported yet and will be skipped.');
end
options_mom_.hessian.use_penalized_objective = false; % penalized objective not yet supported
if isfield(options_,'trend_coeffs')
    error('method_of_moments: %s does not allow for trend in data',mom_method);
end
if ~isempty(options_.endogenous_prior_restrictions.irf) && ~isempty(options_.endogenous_prior_restrictions.moment)
    warning('method_of_moments: Endogenous prior restrictions are not supported yet and will be skipped.');
end
options_mom_.endogenous_prior_restrictions.irf    = {};
options_mom_.endogenous_prior_restrictions.moment = {};
if isfield(options_mom_,'bayesian_irf') && options_mom_.bayesian_irf % do we need this at all??
    warning('method_of_moments: The ''bayesian_irf'' option is not supported yet and will be skipped.');
end
options_mom_.bayesian_irf = false;
if strcmp(mom_method,'IRF_MATCHING')
    if isfield(options_mom_.mom,'penalized_estimator') && options_mom_.mom.penalized_estimator
        warning('method_of_moments: The ''penalized_estimator'' option is not supported yet for IRF_MATCHING and will be ignored.');
    end
    options_mom_.mom.penalized_estimator = false;
end

% -------------------------------------------------------------------------
% OPTIONS POSSIBLY SET BY THE USER
% -------------------------------------------------------------------------

% common settings
options_mom_ = set_default_option(options_mom_,'dirname',dname);         % specify directory in which to store estimation output
options_mom_ = set_default_option(options_mom_,'graph_format','eps');    % specify the file format(s) for graphs saved to disk
options_mom_ = set_default_option(options_mom_,'nodisplay',false);       % do not display the graphs, but still save them to disk
options_mom_ = set_default_option(options_mom_,'nograph',false);         % do not create graphs (which implies that they are not saved to the disk nor displayed)
options_mom_ = set_default_option(options_mom_,'noprint',false);         % do not print output to console
options_mom_ = set_default_option(options_mom_,'TeX',false);             % print TeX tables and graphics
options_mom_.mom = set_default_option(options_mom_.mom,'verbose',false); % display and store intermediate estimation results
if do_bayesian_estimation
    options_mom_ = set_default_option(options_mom_,'plot_priors',true);  % control plotting of priors
    options_mom_ = set_default_option(options_mom_,'prior_trunc',1e-10); % probability of extreme values of the prior density that is ignored when computing bounds for the parameters
end

% specific method_of_moments settings
if strcmp(mom_method,'GMM') || strcmp(mom_method,'SMM')
    options_mom_.mom = set_default_option(options_mom_.mom,'bartlett_kernel_lag',20);                   % bandwith in optimal weighting matrix
    options_mom_.mom = set_default_option(options_mom_.mom,'penalized_estimator',false);                % include deviation from prior mean as additional moment restriction and use prior precision as weights
    options_mom_.mom = set_default_option(options_mom_.mom,'se_tolx',1e-5);                             % step size for numerical computation of standard errors
    options_mom_.mom = set_default_option(options_mom_.mom,'weighting_matrix_scaling_factor',1);        % scaling of weighting matrix in objective function
    options_mom_.mom = set_default_option(options_mom_.mom,'weighting_matrix',{'DIAGONAL'; 'OPTIMAL'}); % weighting matrix in moments distance objective function at each iteration of estimation;
                                                                                                        % possible values are 'OPTIMAL', 'IDENTITY_MATRIX' ,'DIAGONAL' or a filename. Size of cell determines stages in iterated estimation.    
end
if strcmp(mom_method,'SMM')
    options_mom_.mom = set_default_option(options_mom_.mom,'burnin',500);                  % number of periods dropped at beginning of simulation
    options_mom_.mom = set_default_option(options_mom_.mom,'bounded_shock_support',false); % trim shocks in simulation to +- 2 stdev
    options_mom_.mom = set_default_option(options_mom_.mom,'seed',24051986);               % seed used in simulations
    options_mom_.mom = set_default_option(options_mom_.mom,'simulation_multiple',7);       % multiple of the data length used for simulation    
end
if strcmp(mom_method,'GMM')
    options_mom_.mom = set_default_option(options_mom_.mom,'analytic_standard_errors',false); % compute standard errors numerically (0) or analytically (1). Analytical derivatives are only available for GMM.
end
if strcmp(mom_method,'IRF_MATCHING')
    if ~isfield(options_mom_.mom,'irf_matching_file') 
        options_mom_.mom.irf_matching_file = [];  % irf_matching file enables to transform model IRFs before matching them to data IRFs
    end
    options_mom_.mom.irf_matching_file = set_default_option(options_mom_.mom.irf_matching_file,'name','');
    options_mom_.mom = set_default_option(options_mom_.mom,'simulation_method','STOCH_SIMUL'); % simulation method used to compute IRFs
    options_mom_ = set_default_option(options_mom_,'add_tiny_number_to_cholesky',1e-14);       % add tiny number to Cholesky factor to avoid numerical problems when computing IRFs
    options_mom_ = set_default_option(options_mom_,'drop',100);                                % truncation / burnin for order>1 irf simulations
    options_mom_ = set_default_option(options_mom_,'relative_irf',false);                      % requests the computation of normalized IRFs    
end   

% data related options
if strcmp(mom_method,'GMM') || strcmp(mom_method,'SMM')
    options_mom_ = set_default_option(options_mom_,'first_obs',1);     % number of first observation
    options_mom_ = set_default_option(options_mom_,'logdata',false);   % if data is already in logs
    options_mom_ = set_default_option(options_mom_,'nobs',NaN);        % number of observations
    options_mom_ = set_default_option(options_mom_,'prefilter',false); % demean each data series by its empirical mean and use centered moments
    options_mom_ = set_default_option(options_mom_,'xls_sheet',1);     % name of sheet with data in Excel, Octave does not support the empty string, rather use first sheet
    options_mom_ = set_default_option(options_mom_,'xls_range','');    % range of data in Excel sheet    
end

% optimization related
if strcmp(mom_method,'GMM') || strcmp(mom_method,'SMM')
    if (isoctave && user_has_octave_forge_package('optim')) || (~isoctave && user_has_matlab_license('optimization_toolbox'))
        options_mom_ = set_default_option(options_mom_,'mode_compute',13); % specifies lsqnonlin as default optimizer for minimization
    else
        options_mom_ = set_default_option(options_mom_,'mode_compute',5); % specifies newrat as fallback default option for minimization
    end
elseif strcmp(mom_method,'IRF_MATCHING')
    options_mom_ = set_default_option(options_mom_,'mode_compute',5); % specifies newrat as fallback default option for minimization
end
options_mom_ = set_default_option(options_mom_,'additional_optimizer_steps',[]);   % vector of additional mode-finders run after mode_compute
options_mom_ = set_default_option(options_mom_,'optim_opt',[]);                    % a list of NAME and VALUE pairs to set options for the optimization routines. Available options depend on mode_compute
options_mom_ = set_default_option(options_mom_,'silent_optimizer',false);          % run minimization of moments distance silently without displaying results or saving files in between
options_mom_ = set_default_option(options_mom_,'huge_number',1e7);                 % value for replacing the infinite bounds on parameters by finite numbers. Used by some optimizers for numerical reasons
options_mom_.mom = set_default_option(options_mom_.mom,'analytic_jacobian',false); % use analytic Jacobian in optimization, only available for GMM and gradient-based optimizers
options_mom_.optimizer_vec = [options_mom_.mode_compute;num2cell(options_mom_.additional_optimizer_steps)];
options_mom_.mom.analytic_jacobian_optimizers = [1, 3, 4, 13, 101];                % these are currently supported optimizers that are able to use the analytic_jacobian option
options_mom_.analytic_derivation = 0;                                              % force to 0 as we check this seperately in dynare_minimize_objective.m
options_mom_ = set_default_option(options_mom_,'mode_file','');                    % name of the file containing initial values for the mode
options_mom_ = set_default_option(options_mom_,'cova_compute',true);               % 1: computed covariance via Hessian after the computation of the mode, 0: turn off computation of covariance matrix

% perturbation related
options_mom_ = set_default_option(options_mom_,'order',1);                              % order of Taylor approximation in perturbation
if strcmp(mom_method,'IRF_MATCHING')                                                    % number of simulated series used to compute IRFs
    if options_mom_.order == 1
        options_mom_ = set_default_option(options_mom_,'replic',1);
    else
        options_mom_ = set_default_option(options_mom_,'replic',50);
    end
end
options_mom_ = set_default_option(options_mom_,'pruning',false);                        % use pruned state space system at order>1
options_mom_ = set_default_option(options_mom_,'aim_solver',false);                     % use AIM algorithm to compute perturbation approximation instead of mjdgges
options_mom_ = set_default_option(options_mom_,'k_order_solver',false);                 % use k_order_perturbation instead of mjdgges
options_mom_ = set_default_option(options_mom_,'dr_cycle_reduction',false);             % use cycle reduction algorithm to solve the polynomial equation for retrieving the coefficients associated to the endogenous variables in the decision rule
options_mom_ = set_default_option(options_mom_,'dr_cycle_reduction_tol',1e-7);          % convergence criterion used in the cycle reduction algorithm
options_mom_ = set_default_option(options_mom_,'dr_logarithmic_reduction',false);       % use logarithmic reduction algorithm to solve the polynomial equation for retrieving the coefficients associated to the endogenous variables in the decision rule
options_mom_ = set_default_option(options_mom_,'dr_logarithmic_reduction_maxiter',100); % maximum number of iterations used in the logarithmic reduction algorithm
options_mom_ = set_default_option(options_mom_,'dr_logarithmic_reduction_tol',1e-12);   % convergence criterion used in the cycle reduction algorithm
options_mom_ = set_default_option(options_mom_,'qz_criterium',1-1e-6);                  % value used to split stable from unstable eigenvalues in reordering the Generalized Schur decomposition used for solving first order problems
                                                                                        % if there are no unit roots one can use 1.0 (or slightly below) which we set as default; if they are possible, you may have have multiple unit roots and the accuracy decreases when computing the eigenvalues in lyapunov_symm
                                                                                        % Note that unit roots are only possible at first-order, at higher order we set it to 1 in pruned_state_space_system and focus only on stationary observables.
options_mom_ = set_default_option(options_mom_,'qz_zero_threshold',1e-6);               % value used to test if a generalized eigenvalue is 0/0 in the generalized Schur decomposition
options_mom_ = set_default_option(options_mom_,'schur_vec_tol',1e-11);                  % tolerance level used to find nonstationary variables in Schur decomposition of the transition matrix.

% numerical algorithms
options_mom_ = set_default_option(options_mom_,'lyapunov_db',false);                    % doubling algorithm (disclyap_fast) to solve Lyapunov equation to compute variance-covariance matrix of state variables
options_mom_ = set_default_option(options_mom_,'lyapunov_fp',false);                    % fixed-point algorithm to solve Lyapunov equation to compute variance-covariance matrix of state variables
options_mom_ = set_default_option(options_mom_,'lyapunov_srs',false);                   % square-root-solver (dlyapchol) algorithm to solve Lyapunov equation to compute variance-covariance matrix of state variables
options_mom_ = set_default_option(options_mom_,'lyapunov_complex_threshold',1e-15);     % complex block threshold for the upper triangular matrix in symmetric Lyapunov equation solver
options_mom_ = set_default_option(options_mom_,'lyapunov_fixed_point_tol',1e-10);       % convergence criterion used in the fixed point Lyapunov solver
options_mom_ = set_default_option(options_mom_,'lyapunov_doubling_tol',1e-16);          % convergence criterion used in the doubling algorithm

% Bayesian MCMC related
if do_bayesian_estimation
    options_mom_ = set_default_option(options_mom_,'mh_replic',0);                         % number of draws in Metropolis-Hastings and slice samplers
    options_mom_ = set_default_option(options_mom_,'mh_posterior_mode_estimation',false);  % skip optimizer-based mode-finding and instead compute the mode based on a run of a MCMC
    options_mom_ = set_default_option(options_mom_,'load_mh_file',false);                  % add to previous Metropolis-Hastings or slice simulations instead of starting from scratch
    options_mom_ = set_default_option(options_mom_,'load_results_after_load_mh',false);    % load the previously computed convergence diagnostics, marginal data density, and posterior statistics from an existing mom_results file instead of recomputing them
    
    if options_mom_.mh_replic > 0 || options_mom_.load_mh_file        
        options_mom_ = set_default_option(options_mom_,'sub_draws',[]);
        options_mom_ = set_default_option(options_mom_,'posterior_max_subsample_draws',1200);
        options_mom_ = set_default_option(options_mom_,'mh_nblck',2);                         % number of parallel chains for Metropolis-Hastings or slice algorithm
        options_mom_ = set_default_option(options_mom_,'mh_drop',0.5);                        % fraction of initially generated parameter vectors to be dropped as a burn-in before using posterior simulations
        options_mom_ = set_default_option(options_mom_,'mh_conf_sig',0.9);                    % confidence/HPD interval used for the computation of prior and posterior statistics
        options_mom_ = set_default_option(options_mom_,'mh_recover',false);                   % attempts to recover a Metropolis-Hastings simulation that crashed prematurely
        options_mom_ = set_default_option(options_mom_,'MCMC_jumping_covariance','hessian');  % which covariance to use for the proposal density of the MCMC sampler
        if ~isfield(options_mom_,'mh_initialize_from_previous_mcmc')
            options_mom_.mh_initialize_from_previous_mcmc.status = false;                     % pick initial values for new MCMC from a previous one
        end
        options_mom_.mh_initialize_from_previous_mcmc = set_default_option(options_mom_.mh_initialize_from_previous_mcmc,'directory',''); % pick initial values for new MCMC from a previous one: directory
        options_mom_.mh_initialize_from_previous_mcmc = set_default_option(options_mom_.mh_initialize_from_previous_mcmc,'record','');    % pick initial values for new MCMC from a previous one: record file name
        options_mom_.mh_initialize_from_previous_mcmc = set_default_option(options_mom_.mh_initialize_from_previous_mcmc,'prior','');     % pick initial values for new MCMC from a previous one: prior file name
        if ~isfield(options_mom_,'posterior_sampler_options')
            options_mom_.posterior_sampler_options = [];  
        end
        options_mom_.posterior_sampler_options = set_default_option(options_mom_.posterior_sampler_options,'posterior_sampling_method','random_walk_metropolis_hastings'); % selects the sampler used to sample from the posterior distribution during Bayesian estimation
        options_mom_.posterior_sampler_options = set_default_option(options_mom_.posterior_sampler_options,'sampling_opt',[]);                                             % used to set options for the posterior sampling methods
        switch options_mom_.posterior_sampler_options.posterior_sampling_method
            case 'random_walk_metropolis_hastings'
                if ~isfield(options_mom_.posterior_sampler_options,'rwmh')
                    options_mom_.posterior_sampler_options.rwmh = [];
                end
                options_mom_.posterior_sampler_options.rwmh = set_default_option(options_mom_.posterior_sampler_options.rwmh,'proposal_distribution','rand_multivariate_normal');
                options_mom_.posterior_sampler_options.rwmh = set_default_option(options_mom_.posterior_sampler_options.rwmh,'student_degrees_of_freedom',3);
                options_mom_.posterior_sampler_options.rwmh = set_default_option(options_mom_.posterior_sampler_options.rwmh,'use_mh_covariance_matrix',false);
                options_mom_.posterior_sampler_options.rwmh = set_default_option(options_mom_.posterior_sampler_options.rwmh,'save_tmp_file',false);
            case 'tailored_random_block_metropolis_hastings'
                if ~isfield(options_mom_.posterior_sampler_options,'tarb')
                    options_mom_.posterior_sampler_options.tarb = [];
                end
                options_mom_.posterior_sampler_options.tarb = set_default_option(options_mom_.posterior_sampler_options.tarb,'proposal_distribution','rand_multivariate_normal');
                options_mom_.posterior_sampler_options.tarb = set_default_option(options_mom_.posterior_sampler_options.tarb,'student_degrees_of_freedom',3);
                options_mom_.posterior_sampler_options.tarb = set_default_option(options_mom_.posterior_sampler_options.tarb,'mode_compute',4);
                options_mom_.posterior_sampler_options.tarb = set_default_option(options_mom_.posterior_sampler_options.tarb,'new_block_probability',0.25);
                options_mom_.posterior_sampler_options.tarb = set_default_option(options_mom_.posterior_sampler_options.tarb,'optim_opt','');
                options_mom_.posterior_sampler_options.tarb = set_default_option(options_mom_.posterior_sampler_options.tarb,'save_tmp_file',true);
            case 'slice'
                if ~isfield(options_mom_.posterior_sampler_options,'slice')
                    options_mom_.posterior_sampler_options.slice = [];
                end
                options_mom_.posterior_sampler_options.slice = set_default_option(options_mom_.posterior_sampler_options.slice,'proposal_distribution','');
                options_mom_.posterior_sampler_options.slice = set_default_option(options_mom_.posterior_sampler_options.slice,'rotated',0);
                options_mom_.posterior_sampler_options.slice = set_default_option(options_mom_.posterior_sampler_options.slice,'slice_initialize_with_mode',false);  % must be used with rotated
                options_mom_.posterior_sampler_options.slice = set_default_option(options_mom_.posterior_sampler_options.slice,'use_mh_covariance_matrix',false); % must be used with rotated
                options_mom_.posterior_sampler_options.slice = set_default_option(options_mom_.posterior_sampler_options.slice,'WR',[]);        
                options_mom_.posterior_sampler_options.slice = set_default_option(options_mom_.posterior_sampler_options.slice,'mode_files',[]);
                options_mom_.posterior_sampler_options.slice = set_default_option(options_mom_.posterior_sampler_options.slice,'mode',[]);        
                options_mom_.posterior_sampler_options.slice = set_default_option(options_mom_.posterior_sampler_options.slice,'initial_step_size',0.8);
                options_mom_.posterior_sampler_options.slice = set_default_option(options_mom_.posterior_sampler_options.slice,'save_tmp_file',true);
            case 'independent_metropolis_hastings'
                if ~isfield(options_mom_.posterior_sampler_options,'imh')
                    options_mom_.posterior_sampler_options.imh = [];
                end
                options_mom_.posterior_sampler_options.imh = set_default_option(options_mom_.posterior_sampler_options.imh,'proposal_distribution','rand_multivariate_normal');
                options_mom_.posterior_sampler_options.imh = set_default_option(options_mom_.posterior_sampler_options.imh,'use_mh_covariance_matrix',false);
                options_mom_.posterior_sampler_options.imh = set_default_option(options_mom_.posterior_sampler_options.imh,'save_tmp_file',false);                
        end
        if ~strcmp(options_mom_.posterior_sampler_options.posterior_sampling_method,'slice')
            options_mom_ = set_default_option(options_mom_,'mh_init_scale_factor',2);
            options_mom_ = set_default_option(options_mom_,'mh_jscale',[]);
        end
        
        % mh_tune_jscale options
        if strcmp(options_mom_.posterior_sampler_options.posterior_sampling_method,'random_walk_metropolis_hastings')
            if ~isfield(options_mom_,'mh_tune_jscale')
                options_mom_.mh_tune_jscale = [];
            end
            options_mom_.mh_tune_jscale = set_default_option(options_mom_.mh_tune_jscale,'status',false);
            options_mom_.mh_tune_jscale = set_default_option(options_mom_.mh_tune_jscale,'target',0.33);
            options_mom_.mh_tune_jscale = set_default_option(options_mom_.mh_tune_jscale,'guess',[]);
            options_mom_.mh_tune_jscale.maxiter = options_.mh_tune_jscale.maxiter;
            options_mom_.mh_tune_jscale.rho = options_.mh_tune_jscale.rho;
            options_mom_.mh_tune_jscale.stepsize = options_.mh_tune_jscale.stepsize;
            options_mom_.mh_tune_jscale.c1 = options_.mh_tune_jscale.c1;
            options_mom_.mh_tune_jscale.c2 = options_.mh_tune_jscale.c2;
            options_mom_.mh_tune_jscale.c3 = options_.mh_tune_jscale.c3;
        end

        % convergence diagnostics
        options_mom_ = set_default_option(options_mom_,'nodiagnostic',false);
        if ~isfield(options_mom_,'convergence')
            options_mom_.convergence = [];
        end
        if ~isfield(options_mom_.convergence,'geweke')
            options_mom_.convergence.geweke = [];
        end
        if ~isfield(options_mom_.convergence,'rafterylewis')
            options_mom_.convergence.rafterylewis = [];
        end
        if ~isfield(options_mom_.convergence,'brooksgelman')
            options_mom_.convergence.brooksgelman = [];
        end
        options_mom_.convergence.geweke = set_default_option(options_mom_.convergence.geweke,'taper_steps', [4 8 15]);
        options_mom_.convergence.geweke = set_default_option(options_mom_.convergence.geweke,'geweke_interval', [0.2 0.5]);
        options_mom_.convergence.rafterylewis = set_default_option(options_mom_.convergence.rafterylewis,'indicator', false);
        options_mom_.convergence.rafterylewis = set_default_option(options_mom_.convergence.rafterylewis,'qrs', [0.025 0.005 0.95]);
        options_mom_.convergence.brooksgelman = set_default_option(options_mom_.convergence.brooksgelman,'plotrows',3);
    end
end

% mode check plot options
options_mom_.mode_check.nolik = false;                                                          % we don't do likelihood (also this initializes mode_check substructure)
options_mom_.mode_check = set_default_option(options_mom_.mode_check,'status',false);           % plot the target function for values around the computed mode for each estimated parameter in turn. This is helpful to diagnose problems with the optimizer.
options_mom_.mode_check = set_default_option(options_mom_.mode_check,'neighbourhood_size',.5);  % width of the window around the computed mode to be displayed on the diagnostic plots. This width is expressed in percentage deviation. The Inf value is allowed, and will trigger a plot over the entire domain
options_mom_.mode_check = set_default_option(options_mom_.mode_check,'symmetric_plots',true);   % ensure that the check plots are symmetric around the mode. A value of 0 allows to have asymmetric plots, which can be useful if the mode is close to a domain boundary, or in conjunction with neighbourhood_size = Inf when the domain is not the entire real line
options_mom_.mode_check = set_default_option(options_mom_.mode_check,'number_of_points',20);    % number of points around the mode where the target function is evaluated (for each parameter)


% -------------------------------------------------------------------------
% OPTIONS THAT NEED TO BE CARRIED OVER (E.G. SET BY THE PREPROCESSOR)
% -------------------------------------------------------------------------

% related to VAROBS block
options_mom_.varobs = options_.varobs;              % observable variables in order they are declared in varobs
options_mom_.varobs_id = options_.varobs_id;        % index for observable variables in M_.endo_names
options_mom_.obs_nbr = length(options_mom_.varobs); % number of observed variables

% related to call of dynare
options_mom_.console_mode = options_.console_mode;
if options_mom_.console_mode
    options_mom_.nodisplay = true;
end
options_mom_.parallel = options_.parallel;
options_mom_.parallel_info = options_.parallel_info;
options_mom_.debug = options_.debug; % debug option is needed by some functions, e.g. check_plot

% related to estimated_params and estimated_params_init blocks
options_mom_.use_calibration_initialization = options_.use_calibration_initialization;

% related to model block
options_mom_.linear   = options_.linear;
options_mom_.use_dll  = options_.use_dll;
options_mom_.block    = options_.block;
options_mom_.bytecode = options_.bytecode;

% related to steady-state computations
options_mom_.homotopy_force_continue = options_.homotopy_force_continue;
options_mom_.homotopy_mode           = options_.homotopy_mode;
options_mom_.homotopy_steps          = options_.homotopy_steps;
options_mom_.markowitz               = options_.markowitz;
options_mom_.solve_algo              = options_.solve_algo;
options_mom_.solve_tolf              = options_.solve_tolf;
options_mom_.solve_tolx              = options_.solve_tolx;
options_mom_.steady                  = options_.steady;
options_mom_.steadystate             = options_.steadystate;
options_mom_.steadystate_flag        = options_.steadystate_flag;
options_mom_.threads = options_.threads; % needed by resol
options_mom_.debug = options_.debug; % debug option needed by some functions, e.g. check_plot

% random numbers
options_mom_.DynareRandomStreams.seed = options_.DynareRandomStreams.seed;
options_mom_.DynareRandomStreams.algo = options_.DynareRandomStreams.algo;

% dataset_ related
options_mom_.dataset        = options_.dataset;
options_mom_.initial_period = options_.initial_period;

% optimization related
if any(cellfun(@(x) isnumeric(x) && any(x == 2), options_mom_.optimizer_vec)) % simulated annealing (mode_compute=2)
    options_mom_.saopt = options_.saopt;
end
if any(cellfun(@(x) isnumeric(x) && any(x == 4), options_mom_.optimizer_vec)) % csminwel (mode_compute=4)
    options_mom_.csminwel = options_.csminwel;
end
if any(cellfun(@(x) isnumeric(x) && any(x == 5), options_mom_.optimizer_vec)) % newrat (mode_compute=5)
    options_mom_.newrat = options_.newrat;
end
if any(cellfun(@(x) isnumeric(x) && any(x == 6), options_mom_.optimizer_vec)) % gmhmaxlik (mode_compute=6)
    options_mom_.gmhmaxlik = options_.gmhmaxlik;
    options_mom_.mh_jscale = options_.mh_jscale;
end
if any(cellfun(@(x) isnumeric(x) && any(x == 8), options_mom_.optimizer_vec)) % simplex variation on Nelder Mead algorithm (mode_compute=8)
    options_mom_.simplex = options_.simplex;
end
if any(cellfun(@(x) isnumeric(x) && any(x == 9), options_mom_.optimizer_vec)) % cmaes (mode_compute=9)
    options_mom_.cmaes = options_.cmaes;
end
if any(cellfun(@(x) isnumeric(x) && any(x == 10), options_mom_.optimizer_vec)) % simpsa (mode_compute=10)
    options_mom_.simpsa = options_.simpsa;
end
if any(cellfun(@(x) isnumeric(x) && any(x == 12), options_mom_.optimizer_vec)) % particleswarm (mode_compute=12)
    options_mom_.particleswarm = options_.particleswarm;
end
if any(cellfun(@(x) isnumeric(x) && any(x == 101), options_mom_.optimizer_vec)) % solveopt (mode_compute=101)
    options_mom_.solveopt = options_.solveopt;
end
if any(cellfun(@(x) isnumeric(x) && (any(x == 4) || any(x == 5)), options_mom_.optimizer_vec)) % used by csminwel and newrat
    options_mom_.gradient_method = options_.gradient_method;
    options_mom_.gradient_epsilon = options_.gradient_epsilon;
end
options_mom_.gstep = options_.gstep; % needed by hessian.m
options_mom_.trust_region_initial_step_bound_factor = options_.trust_region_initial_step_bound_factor; % used in dynare_solve for trust_region

% miscellaneous
options_mom_.threads = options_.threads;
options_mom_.MaxNumberOfBytes = options_.MaxNumberOfBytes;
options_mom_.marginal_data_density = options_.marginal_data_density;


% -------------------------------------------------------------------------
% DEFAULT VALUES
% -------------------------------------------------------------------------
options_mom_.mom.compute_derivs = false;   % flag to compute derivs in objective function (might change for GMM with either analytic_standard_errors or analytic_jacobian (dependent on optimizer))
options_mom_.mom.vector_output = false;    % specifies whether the objective function returns a vector
options_mom_.analytic_derivation_mode = 0; % needed by get_perturbation_params_derivs.m, ie use efficient sylvester equation method to compute analytical derivatives as in Ratto & Iskrev (2012)
options_mom_.initialize_estimated_parameters_with_the_prior_mode = 0; % needed by set_prior.m
options_mom_.figures = options_.figures;   % needed by plot_priors.m
options_mom_.ramsey_policy = false;        % needed by evaluate_steady_state
options_mom_.risky_steadystate = false;    % needed by resol
options_mom_.jacobian_flag = true;         % needed by dynare_solve
options_mom_.use_mh_covariance_matrix = false; % needed by posterior_sampler, get's overwritten by same option in options_mom_.posterior_sampler_options


% -------------------------------------------------------------------------
% CHECKS ON SETTINGS
% -------------------------------------------------------------------------
if strcmp(mom_method,'GMM') || strcmp(mom_method,'SMM')
    if numel(options_mom_.nobs) > 1
        error('method_of_moments: Recursive estimation is not supported. Please set an integer as ''nobs''!');
    end
    if numel(options_mom_.first_obs) > 1
        error('method_of_moments: Recursive estimation is not supported. Please set an integer as ''first_obs''!');
    end
end
if options_mom_.order < 1
    error('method_of_moments: The order of the Taylor approximation cannot be 0!')
end
if options_mom_.order > 2
    fprintf('Dynare will use ''k_order_solver'' as the order>2\n');
    options_mom_.k_order_solver = true;
end
if strcmp(mom_method,'SMM')
    if options_mom_.mom.simulation_multiple < 1
        fprintf('The simulation horizon is shorter than the data. Dynare resets the simulation_multiple to 7.\n')
        options_mom_.mom.simulation_multiple = 7;
    end
end
if strcmp(mom_method,'GMM')
    % require pruning with GMM at higher order
    if options_mom_.order > 1 && ~options_mom_.pruning
        fprintf('GMM at higher order only works with pruning, so we set pruning option to 1.\n');
        options_mom_.pruning = true;
    end
    if options_mom_.order > 3
        error('method_of_moments: Perturbation orders higher than 3 are not implemented for GMM estimation, try using SMM!');
    end
end
if strcmp(mom_method,'IRF_MATCHING') && do_bayesian_estimation
    if isfield(options_mom_,'mh_tune_jscale') && options_mom_.mh_tune_jscale.status && (options_mom_.mh_tune_jscale.maxiter<options_mom_.mh_tune_jscale.stepsize)
        warning('method_of_moments: You specified mh_tune_jscale, but the maximum number of iterations is smaller than the step size. No update will take place.')
    end
    if options_mom_.load_results_after_load_mh
        if ~exist([options_mom_.dirname filesep 'method_of_moments' filesep M_.fname '_mom_results.mat'],'file')
            fprintf('\nYou specified the ''load_results_after_load_mh'' option, but no ''%s_mom_results.mat'' file\n',M_.fname);
            fprintf('was found in the folder %s%smethod_of_moments.\n',options_mom_.dirname,filesep);
            fprintf('Results will be recomputed and option ''load_results_after_load_mh'' is reset to false.\n');
            options_mom_.load_results_after_load_mh = false;
        end
    end
    if options_mom_.mh_replic>0 && options_mom_.mh_nblck<1
        error('method_of_moments: Bayesian MCMC estimation cannot be conducted with ''mh_nblocks''=0!')
    end
end
if options_mom_.mom.analytic_jacobian && ~strcmp(mom_method,'GMM')
    options_mom_.mom.analytic_jacobian = false;
    fprintf('\n''analytic_jacobian'' option will be dismissed as it only works with GMM.\n');
end
if strcmp(options_mom_.mom.mom_method,'IRF_MATCHING')
    if any(cellfun(@(x) isnumeric(x) && any(x == 13), options_mom_.optimizer_vec))
        error('method_of_moments: lsqnonlin (mode_compute=13) is not yet supported for IRF Matching!');
    end
end