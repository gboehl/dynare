function [oo_, options_mom_, M_] = run(bayestopt_, options_, oo_, estim_params_, M_, options_mom_)
%function [oo_, options_mom_, M_] = run(bayestopt_, options_, oo_, estim_params_, M_, options_mom_)
% -------------------------------------------------------------------------
% This function performs a method of moments estimation with the following steps:
%   Step 0: Check if required structures and options exist
%   Step 1: - Prepare options_mom_ structure
%           - Carry over options from the preprocessor
%           - Initialize other options
%           - Get variable orderings and state space representation
%   Step 2: Checks and transformations for matched moments structure
%   Step 3: Checks and transformations for estimated parameters, priors, and bounds
%   Step 4: Checks and transformations for data
%   Step 5: Checks for steady state at initial parameters
%   Step 6: Checks for objective function at initial parameters
%   Step 7: Iterated method of moments estimation
%   Step 8: J-Test
%   Step 9: Clean up
% -------------------------------------------------------------------------
% This function is inspired by replication codes accompanied to the following papers:
%  o Andreasen, Fernández-Villaverde, Rubio-Ramírez (2018): "The Pruned State-Space System for Non-Linear DSGE Models: Theory and Empirical Applications", Review of Economic Studies, 85(1):1-49.
%  o Born, Pfeifer (2014): "Risk Matters: Comment", American Economic Review, 104(12):4231-4239.
%  o Mutschler (2018): "Higher-order statistics for DSGE models", Econometrics and Statistics, 6:44-56.
% =========================================================================
% INPUTS
%  o bayestopt_:             [structure] information about priors
%  o options_:               [structure] information about global options
%  o oo_:                    [structure] storage for results
%  o estim_params_:          [structure] information about estimated parameters
%  o M_:                     [structure] information about model with
%                               o matched_moments:       [cell] information about selected moments to match in estimation
%                                                               vars: matched_moments{:,1});
%                                                               lead/lags: matched_moments{:,2}; 
%                                                               powers: matched_moments{:,3};
%  o options_mom_:           [structure] information about settings specified by the user
% -------------------------------------------------------------------------
% OUTPUTS
%  o oo_:                    [structure] storage for results (oo_)
%  o options_mom_:           [structure] information about all (user-specified and updated) settings used in estimation (options_mom_)
% -------------------------------------------------------------------------
% This function is called by
%  o driver.m
% -------------------------------------------------------------------------
% This function calls
%  o check_for_calibrated_covariances.m
%  o check_prior_bounds.m
%  o do_parameter_initialization.m
%  o dynare_minimize_objective.m
%  o evaluate_steady_state
%  o get_all_parameters.m
%  o get_matrix_entries_for_psd_check.m
%  o makedataset.m
%  o mom.check_plot.m
%  o mom.data_moments.m
%  o mom.objective_function.m
%  o mom.optimal_weighting_matrix
%  o mom-standard_errors
%  o plot_priors.m
%  o print_info.m
%  o prior_bounds.m
%  o set_default_option.m
%  o set_prior.m
%  o set_state_space.m
%  o set_all_parameters.m
%  o test_for_deep_parameters_calibration.m
% =========================================================================
% Copyright (C) 2020-2022 Dynare Team
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
% Author(s): 
% o Willi Mutschler (willi@mutschler.eu)
% o Johannes Pfeifer (jpfeifer@uni-koeln.de)
% =========================================================================

%% TO DO LIST
% - [ ] add IRF matching
% - [ ] speed up pruned_state_space_system (by using doubling with old initial values, hardcoding zeros, other "tricks" used in e.g. nlma)
% - [ ] add option to use autocorrelations (we have useautocorr in identification toolbox already)
% - [ ] SMM with extended path
% - [ ] deal with measurement errors (once @wmutschl has implemented this in identification toolbox)
% - [ ] dirname option to save output to different directory not yet implemented
% - [ ] display scaled moments

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

% -------------------------------------------------------------------------
% Step 0: Check if required structures and options exist
% -------------------------------------------------------------------------
if isempty(estim_params_) % structure storing the info about estimated parameters in the estimated_params block
    if ~(isfield(estim_params_,'nvx') && (size(estim_params_.var_exo,1)+size(estim_params_.var_endo,1)+size(estim_params_.corrx,1)+size(estim_params_.corrn,1)+size(estim_params_.param_vals,1))==0)
        error('method_of_moments: You need to provide an ''estimated_params'' block')
    else
        error('method_of_moments: The ''estimated_params'' block must not be empty')
    end
end
if ~isfield(M_,'matched_moments') || isempty(M_.matched_moments) % structure storing the moments used for the method of moments estimation
    error('method_of_moments: You need to provide a ''matched_moments'' block')
end
if ~isempty(bayestopt_) && any(bayestopt_.pshape==0) && any(bayestopt_.pshape~=0)
    error('method_of_moments: Estimation must be either fully classical or fully Bayesian. Maybe you forgot to specify a prior distribution.')
end

if options_.logged_steady_state || options_.loglinear
    error('method_of_moments: The loglinear option is not supported. Please append the required logged variables as auxiliary equations.\n')
else
    options_mom_.logged_steady_state = 0;
    options_mom_.loglinear = false;
end

fprintf('\n==== Method of Moments Estimation (%s) ====\n\n',options_mom_.mom.mom_method)

% -------------------------------------------------------------------------
% Step 1a: Prepare options_mom_ structure
% -------------------------------------------------------------------------
% options_mom_ is local and contains default and user-specified values for 
% all settings needed for the method of moments estimation. Some options,
% though, are set by the preprocessor into options_ and we copy these over.
% The idea is to be independent of options_ and have full control of the
% estimation instead of possibly having to deal with options chosen somewhere
% else in the mod file.

% Method of Moments estimation options that can be set by the user in the mod file, otherwise default values are provided
if strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')
    options_mom_.mom = set_default_option(options_mom_.mom,'bartlett_kernel_lag',20);               % bandwith in optimal weighting matrix
    options_mom_.mom = set_default_option(options_mom_.mom,'penalized_estimator',false);            % include deviation from prior mean as additional moment restriction and use prior precision as weight
    options_mom_.mom = set_default_option(options_mom_.mom,'verbose',false);                        % display and store intermediate estimation results
    options_mom_.mom = set_default_option(options_mom_.mom,'weighting_matrix',{'DIAGONAL'; 'OPTIMAL'});   % weighting matrix in moments distance objective function at each iteration of estimation;
                                                                                                           % possible values are 'OPTIMAL', 'IDENTITY_MATRIX' ,'DIAGONAL' or a filename. Size of cell determines stages in iterated estimation.
    options_mom_.mom = set_default_option(options_mom_.mom,'weighting_matrix_scaling_factor',1);    % scaling of weighting matrix in objective function
    options_mom_.mom = set_default_option(options_mom_.mom,'se_tolx',1e-5);                         % step size for numerical computation of standard errors
    options_mom_ = set_default_option(options_mom_,'order',1);                                      % order of Taylor approximation in perturbation
    options_mom_ = set_default_option(options_mom_,'pruning',false);                                % use pruned state space system at higher-order
    % Checks for perturbation order
    if options_mom_.order < 1
        error('method_of_moments: The order of the Taylor approximation cannot be 0!')
    end
end
if strcmp(options_mom_.mom.mom_method,'SMM')
    options_mom_.mom = set_default_option(options_mom_.mom,'burnin',500);                           % number of periods dropped at beginning of simulation
    options_mom_.mom = set_default_option(options_mom_.mom,'bounded_shock_support',false);          % trim shocks in simulation to +- 2 stdev
    options_mom_.mom = set_default_option(options_mom_.mom,'seed',24051986);                        % seed used in simulations
    options_mom_.mom = set_default_option(options_mom_.mom,'simulation_multiple',7);                % multiple of the data length used for simulation
    if options_mom_.mom.simulation_multiple < 1
        fprintf('The simulation horizon is shorter than the data. Dynare resets the simulation_multiple to 5.\n')
        options_mom_.mom.simulation_multiple = 7;
    end
end
if strcmp(options_mom_.mom.mom_method,'GMM')
    % Check for pruning with GMM at higher order
    if options_mom_.order > 1 && ~options_mom_.pruning
        fprintf('GMM at higher order only works with pruning, so we set pruning option to 1.\n');
        options_mom_.pruning = true;
    end
    if options_mom_.order > 3
        error('method_of_moments: perturbation orders higher than 3 are not implemented for GMM estimation, try using SMM.\n');
    end
end
options_mom_.mom = set_default_option(options_mom_.mom,'analytic_standard_errors',false);       % compute standard errors numerically (0) or analytically (1). Analytical derivatives are only available for GMM.
options_mom_.mom = set_default_option(options_mom_.mom,'analytic_jacobian',false);              % use analytic Jacobian in optimization, only available for GMM and gradient-based optimizers
% initialize flag to compute derivs in objective function (needed for GMM with either analytic_standard_errors or analytic_jacobian )
options_mom_.mom.compute_derivs = false;
    
% General options that can be set by the user in the mod file, otherwise default values are provided
options_mom_ = set_default_option(options_mom_,'dirname',M_.dname);    % specify directory in which to store estimation output [not yet working]
options_mom_ = set_default_option(options_mom_,'graph_format','eps');  % specify the file format(s) for graphs saved to disk
options_mom_ = set_default_option(options_mom_,'nodisplay',false);     % do not display the graphs, but still save them to disk
options_mom_ = set_default_option(options_mom_,'nograph',false);       % do not create graphs (which implies that they are not saved to the disk nor displayed)
options_mom_ = set_default_option(options_mom_,'noprint',false);       % do not print output to console
options_mom_ = set_default_option(options_mom_,'plot_priors',true);    % control plotting of priors
options_mom_ = set_default_option(options_mom_,'prior_trunc',1e-10);   % probability of extreme values of the prior density that is ignored when computing bounds for the parameters
options_mom_ = set_default_option(options_mom_,'TeX',false);           % print TeX tables and graphics
options_mom_ = set_default_option(options_mom_,'verbosity',false);           % print TeX tables and graphics

% Data and model options that can be set by the user in the mod file, otherwise default values are provided
options_mom_ = set_default_option(options_mom_,'first_obs',1);     % number of first observation
options_mom_ = set_default_option(options_mom_,'logdata',false);   % if data is already in logs
options_mom_ = set_default_option(options_mom_,'nobs',NaN);        % number of observations
options_mom_ = set_default_option(options_mom_,'prefilter',false); % demean each data series by its empirical mean and use centered moments
options_mom_ = set_default_option(options_mom_,'xls_sheet',1);     % name of sheet with data in Excel
options_mom_ = set_default_option(options_mom_,'xls_range','');    % range of data in Excel sheet
% temporary workaround for https://git.dynare.org/Dynare/dseries/-/issues/51
if options_mom_.xls_sheet~=1
    evalin('base','options_.xls_sheet=options_mom_.xls_sheet');
end
if ~isempty(options_mom_.xls_range)
    evalin('base','options_.xls_range=options_mom_.xls_range');
end

% Recursive estimation and forecast are not supported
if numel(options_mom_.nobs)>1
    error('method_of_moments: Recursive estimation and forecast for samples is not supported. Please set an integer as ''nobs''.');
end
if numel(options_mom_.first_obs)>1
    error('method_of_moments: Recursive estimation and forecast for samples is not supported. Please set an integer as ''first_obs''.');
end

% Optimization options that can be set by the user in the mod file, otherwise default values are provided
options_mom_ = set_default_option(options_mom_,'huge_number',1e7);               % value for replacing the infinite bounds on parameters by finite numbers. Used by some optimizers for numerical reasons
if (isoctave && user_has_octave_forge_package('optim')) || (~isoctave && user_has_matlab_license('optimization_toolbox'))
    options_mom_ = set_default_option(options_mom_,'mode_compute',13);               % specifies lsqnonlin as default optimizer for minimization of moments distance
else
    options_mom_ = set_default_option(options_mom_,'mode_compute',4);               % specifies csminwel as fallback default option for minimization of moments distance
end
options_mom_ = set_default_option(options_mom_,'additional_optimizer_steps',[]); % vector of additional mode-finders run after mode_compute
options_mom_ = set_default_option(options_mom_,'optim_opt',[]);                  % a list of NAME and VALUE pairs to set options for the optimization routines. Available options depend on mode_compute
options_mom_ = set_default_option(options_mom_,'silent_optimizer',false);        % run minimization of moments distance silently without displaying results or saving files in between
% Check plot options that can be set by the user in the mod file, otherwise default values are provided
options_mom_.mode_check.nolik = false;                                                          % we don't do likelihood (also this initializes mode_check substructure)
options_mom_.mode_check = set_default_option(options_mom_.mode_check,'status',false);           % plot the target function for values around the computed minimum for each estimated parameter in turn. This is helpful to diagnose problems with the optimizer.
options_mom_.mode_check = set_default_option(options_mom_.mode_check,'neighbourhood_size',.5);  % width of the window around the computed minimum to be displayed on the diagnostic plots. This width is expressed in percentage deviation. The Inf value is allowed, and will trigger a plot over the entire domain
options_mom_.mode_check = set_default_option(options_mom_.mode_check,'symmetric_plots',true);   % ensure that the check plots are symmetric around the minimum. A value of 0 allows to have asymmetric plots, which can be useful if the minimum is close to a domain boundary, or in conjunction with neighbourhood_size = Inf when the domain is not the entire real line
options_mom_.mode_check = set_default_option(options_mom_.mode_check,'number_of_points',20);    % number of points around the minimum where the target function is evaluated (for each parameter)

% Numerical algorithms options that can be set by the user in the mod file, otherwise default values are provided
options_mom_ = set_default_option(options_mom_,'aim_solver',false);                     % use AIM algorithm to compute perturbation approximation instead of mjdgges
options_mom_ = set_default_option(options_mom_,'k_order_solver',false);                 % use k_order_perturbation instead of mjdgges
options_mom_ = set_default_option(options_mom_,'dr_cycle_reduction',false);             % use cycle reduction algorithm to solve the polynomial equation for retrieving the coefficients associated to the endogenous variables in the decision rule
options_mom_ = set_default_option(options_mom_,'dr_cycle_reduction_tol',1e-7);          % convergence criterion used in the cycle reduction algorithm
options_mom_ = set_default_option(options_mom_,'dr_logarithmic_reduction',false);       % use logarithmic reduction algorithm to solve the polynomial equation for retrieving the coefficients associated to the endogenous variables in the decision rule
options_mom_ = set_default_option(options_mom_,'dr_logarithmic_reduction_maxiter',100); % maximum number of iterations used in the logarithmic reduction algorithm
options_mom_ = set_default_option(options_mom_,'dr_logarithmic_reduction_tol',1e-12);   % convergence criterion used in the cycle reduction algorithm
options_mom_ = set_default_option(options_mom_,'lyapunov_db',false);                    % doubling algorithm (disclyap_fast) to solve Lyapunov equation to compute variance-covariance matrix of state variables
options_mom_ = set_default_option(options_mom_,'lyapunov_fp',false);                    % fixed-point algorithm to solve Lyapunov equation to compute variance-covariance matrix of state variables
options_mom_ = set_default_option(options_mom_,'lyapunov_srs',false);                   % square-root-solver (dlyapchol) algorithm to solve Lyapunov equation to compute variance-covariance matrix of state variables
options_mom_ = set_default_option(options_mom_,'lyapunov_complex_threshold',1e-15);     % complex block threshold for the upper triangular matrix in symmetric Lyapunov equation solver
options_mom_ = set_default_option(options_mom_,'lyapunov_fixed_point_tol',1e-10);       % convergence criterion used in the fixed point Lyapunov solver
options_mom_ = set_default_option(options_mom_,'lyapunov_doubling_tol',1e-16);          % convergence criterion used in the doubling algorithm
options_mom_ = set_default_option(options_mom_,'sylvester_fp',false);                   % determines whether to use fixed point algorihtm to solve Sylvester equation (gensylv_fp), faster for large scale models
options_mom_ = set_default_option(options_mom_,'sylvester_fixed_point_tol',1e-12);      % convergence criterion used in the fixed point Sylvester solver
options_mom_ = set_default_option(options_mom_,'qz_criterium',1-1e-6);                  % value used to split stable from unstable eigenvalues in reordering the Generalized Schur decomposition used for solving first order problems
                                                                                        % if there are no unit roots one can use 1.0 (or slightly below) which we set as default; if they are possible, you may have have multiple unit roots and the accuracy decreases when computing the eigenvalues in lyapunov_symm
                                                                                        % Note that unit roots are only possible at first-order, at higher order we set it to 1 in pruned_state_space_system and focus only on stationary observables.
options_mom_ = set_default_option(options_mom_,'qz_zero_threshold',1e-6);               % value used to test if a generalized eigenvalue is 0/0 in the generalized Schur decomposition
options_mom_ = set_default_option(options_mom_,'schur_vec_tol',1e-11);                  % tolerance level used to find nonstationary variables in Schur decomposition of the transition matrix.
if options_mom_.order > 2
    fprintf('Dynare will use ''k_order_solver'' as the order>2\n');
    options_mom_.k_order_solver = true;
end

% -------------------------------------------------------------------------
% Step 1b: Options that are set by the preprocessor and need to be carried over
% -------------------------------------------------------------------------

% options related to VAROBS
if ~isfield(options_,'varobs')
    error('method_of_moments: VAROBS statement is missing!')
else
    options_mom_.varobs  = options_.varobs;             % observable variables in declaration order
    options_mom_.obs_nbr = length(options_mom_.varobs); % number of observed variables
    % Check that each declared observed variable is also an endogenous variable
    for i = 1:options_mom_.obs_nbr
        if ~any(strcmp(options_mom_.varobs{i},M_.endo_names))
            error(['method_of_moments: Unknown variable (' options_mom_.varobs{i} ')!'])
        end
    end

    % Check that a variable is not declared as observed more than once
    if length(unique(options_mom_.varobs))<length(options_mom_.varobs)
        for i = 1:options_mom_.obs_nbr
            if sum(strcmp(options_mom_.varobs{i},options_mom_.varobs))>1
                error(['method_of_moments: A variable cannot be declared as observed more than once (' options_mom_.varobs{i} ')!'])
            end
        end
    end
end

% options related to variable declarations
if isfield(options_,'trend_coeffs')
    error('method_of_moments: %s does not allow for trend in data',options_mom_.mom.mom_method)
end

% options related to estimated_params and estimated_params_init
options_mom_.use_calibration_initialization = options_.use_calibration_initialization;

% options related to model block
options_mom_.linear   = options_.linear;
options_mom_.use_dll  = options_.use_dll;
options_mom_.block    = options_.block;
options_mom_.bytecode = options_.bytecode;

% options related to steady command
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

% options related to dataset
options_mom_.dataset        = options_.dataset;
options_mom_.initial_period = options_.initial_period;

% options related to endogenous prior restrictions are not supported
options_mom_.endogenous_prior_restrictions.irf    = {};
options_mom_.endogenous_prior_restrictions.moment = {};
if ~isempty(options_.endogenous_prior_restrictions.irf) && ~isempty(options_.endogenous_prior_restrictions.moment)
    fprintf('Endogenous prior restrictions are not supported yet and will be skipped.\n')
end

% -------------------------------------------------------------------------
% Step 1c: Options related to optimizers
% -------------------------------------------------------------------------
% mode_compute = 1, 3, 7, 11, 102, 11, 13
% nothing to be done
% mode_compute = 2
options_mom_.saopt            = options_.saopt;
% mode_compute = 4
options_mom_.csminwel         = options_.csminwel;
% mode_compute = 5
options_mom_.newrat           = options_.newrat;
options_mom_.gstep            = options_.gstep;
% mode_compute = 6
options_mom_.gmhmaxlik        = options_.gmhmaxlik;
options_mom_.mh_jscale        = options_.mh_jscale;
% mode_compute = 8
options_mom_.simplex          = options_.simplex;
% mode_compute = 9
options_mom_.cmaes            = options_.cmaes;
% mode_compute = 10
options_mom_.simpsa           = options_.simpsa;
% mode_compute = 12
options_mom_.particleswarm    = options_.particleswarm;
% mode_compute = 101
options_mom_.solveopt         = options_.solveopt;

options_mom_.gradient_method  = options_.gradient_method;
options_mom_.gradient_epsilon = options_.gradient_epsilon;
options_mom_.analytic_derivation = 0;
options_mom_.analytic_derivation_mode = 0; % needed by get_perturbation_params_derivs.m, ie use efficient sylvester equation method to compute analytical derivatives as in Ratto & Iskrev (2012)

options_mom_.vector_output= false;           % specifies whether the objective function returns a vector

optimizer_vec=[options_mom_.mode_compute;num2cell(options_mom_.additional_optimizer_steps)]; % at each stage one can possibly use different optimizers sequentially

analytic_jacobian_optimizers = [1, 3, 4, 13, 101]; %these are currently supported, see to-do list

% -------------------------------------------------------------------------
% Step 1d: Other options that need to be initialized
% -------------------------------------------------------------------------
options_mom_.initialize_estimated_parameters_with_the_prior_mode = 0; % needed by set_prior.m
options_mom_.figures.textwidth = 0.8; %needed by plot_priors.m
options_mom_.ramsey_policy = 0; % needed by evaluate_steady_state
options_mom_ = set_default_option(options_mom_,'debug',false); %neeeded by e.g. check_plot
options_mom_.risky_steadystate = false; %needed by resol
options_mom_.threads = options_.threads; %needed by resol
options_mom_.jacobian_flag = true;
options_mom_.gstep = options_.gstep;

% options_mom.dsge_var          = false; %needed by check_list_of_variables
% options_mom.bayesian_irf      = false; %needed by check_list_of_variables
% options_mom.moments_varendo   = false; %needed by check_list_of_variables
% options_mom.smoother          = false; %needed by check_list_of_variables
% options_mom.filter_step_ahead = [];  %needed by check_list_of_variables
% options_mom.forecast = 0;
%options_mom_ = set_default_option(options_mom_,'endo_vars_for_moment_computations_in_estimation',[]);

% -------------------------------------------------------------------------
% Step 1e: Get variable orderings and state space representation
% -------------------------------------------------------------------------
oo_.dr = set_state_space(oo_.dr,M_,options_mom_);
% Get index of observed variables in DR order
oo_.dr.obs_var = [];
for i=1:options_mom_.obs_nbr
    oo_.dr.obs_var = [oo_.dr.obs_var; find(strcmp(options_mom_.varobs{i}, M_.endo_names(oo_.dr.order_var)))];
end

% -------------------------------------------------------------------------
% Step 2: Checks and transformations for matched_moments structure
% -------------------------------------------------------------------------
M_.matched_moments_orig = M_.matched_moments; %save original structure

% higher-order product moments not supported yet for GMM
if strcmp(options_mom_.mom.mom_method, 'GMM') && any(cellfun(@sum,M_.matched_moments(:,3))> 2)
    error('method_of_moments: GMM does not yet support product moments higher than 2. Change row %d in ''matched_moments'' block.',jm);
end

% check for duplicate moment conditions
for jm=1:size(M_.matched_moments,1)
    % expand powers to vector of ones
    if any(M_.matched_moments{jm,3}>1)
        tmp1=[]; tmp2=[]; tmp3=[];
        for jjm=1:length(M_.matched_moments{jm,3})
            tmp1 = [tmp1 repmat(M_.matched_moments{jm,1}(jjm),[1 M_.matched_moments{jm,3}(jjm)]) ];
            tmp2 = [tmp2 repmat(M_.matched_moments{jm,2}(jjm),[1 M_.matched_moments{jm,3}(jjm)]) ];
            tmp3 = [tmp3 repmat(1,[1 M_.matched_moments{jm,3}(jjm)]) ];
        end
        M_.matched_moments{jm,1} = tmp1;
        M_.matched_moments{jm,2} = tmp2;
        M_.matched_moments{jm,3} = tmp3;
    end
    % shift time structure to focus only on lags
    M_.matched_moments{jm,2} = M_.matched_moments{jm,2} - max(M_.matched_moments{jm,2});
    % sort such that t=0 variable comes first
    [M_.matched_moments{jm,2},idx_sort] = sort(M_.matched_moments{jm,2},'descend');
    M_.matched_moments{jm,1} = M_.matched_moments{jm,1}(idx_sort);
    M_.matched_moments{jm,3} = M_.matched_moments{jm,3}(idx_sort);
end

% find duplicate rows in cell array by making groups according to powers as we can then use cell2mat for the unique function
powers = cellfun(@sum,M_.matched_moments(:,3))';
UniqueMomIdx = [];
for jpow = unique(powers)
    idx1 = find(powers==jpow);
    [~,idx2] = unique(cell2mat(M_.matched_moments(idx1,:)),'rows');
    UniqueMomIdx = [UniqueMomIdx idx1(idx2)];
end

% remove duplicate elements
DuplicateMoms = setdiff(1:size(M_.matched_moments_orig,1),UniqueMomIdx);
if ~isempty(DuplicateMoms)
    fprintf('Found and removed duplicate declared moments in ''matched_moments'' block in rows:\n %s.\n',num2str(DuplicateMoms))
    fprintf('Dynare will continue with remaining moment conditions\n');
end

if strcmp(options_mom_.mom.mom_method, 'SMM')
    % for SMM we can keep the original structure but get rid of duplicate moments
    M_.matched_moments = M_.matched_moments_orig(sort(UniqueMomIdx),:);
elseif strcmp(options_mom_.mom.mom_method, 'GMM')
    % for GMM we use the transformed matched_moments structure
    M_.matched_moments = M_.matched_moments(sort(UniqueMomIdx),:);
end

% Check if both prefilter and first moments were specified
first_moment_indicator = find(cellfun(@(x) sum(abs(x))==1,M_.matched_moments(:,3)))';
if options_mom_.prefilter && ~isempty(first_moment_indicator)
    fprintf('Centered moments requested (prefilter option is set); therefore, ignore declared first moments in ''matched_moments'' block.\n');
    M_.matched_moments(first_moment_indicator,:)=[]; %remove first moments entries
end
options_mom_.mom.mom_nbr = size(M_.matched_moments,1);

% Get maximum lag number for autocovariances/autocorrelations
options_mom_.ar = max(cellfun(@max,M_.matched_moments(:,2))) - min(cellfun(@min,M_.matched_moments(:,2)));

%check that only observed variables are involved in moments
not_observed_variables=setdiff(oo_.dr.inv_order_var([M_.matched_moments{:,1}]),oo_.dr.obs_var);
if ~isempty(not_observed_variables)
    error('\nmethod_of_moments: You specified moments involving %s, but it is not a varobs.',M_.endo_names{oo_.dr.order_var(not_observed_variables)})
end

% -------------------------------------------------------------------------
% Step 3: Checks and transformations for estimated parameters, priors, and bounds
% -------------------------------------------------------------------------

% Set priors and bounds over the estimated parameters
[xparam0, estim_params_, bayestopt_, lb, ub, M_] = set_prior(estim_params_, M_, options_mom_);

% Check measurement errors
if (estim_params_.nvn || estim_params_.ncn) && strcmp(options_mom_.mom.mom_method, 'GMM')
    error('method_of_moments: GMM estimation does not support measurement error(s) yet. Please specifiy them as a structural shock.')
end

% Check if enough moments for estimation
if options_mom_.mom.mom_nbr < length(xparam0)
    fprintf('\n');
    error('method_of_moments: We must have at least as many moments as parameters for a method of moments estimation.')
end
fprintf('\n\n')

% Check if a _prior_restrictions.m file exists
if exist([M_.fname '_prior_restrictions.m'],'file')
    options_mom_.prior_restrictions.status = 1;
    options_mom_.prior_restrictions.routine = str2func([M_.fname '_prior_restrictions']);
end

bayestopt_laplace=bayestopt_;

% Check on specified priors and penalized estimation
if any(bayestopt_laplace.pshape > 0) % prior specified
    if ~options_mom_.mom.penalized_estimator
        fprintf('\nPriors were specified, but the penalized_estimator-option was not set.\n')
        fprintf('Dynare sets penalized_estimator to 1. Conducting %s with penalty.\n',options_mom_.mom.mom_method)
        options_mom_.mom.penalized_estimator=1;
    end
    if any(setdiff([0;bayestopt_laplace.pshape],[0,3]))
        fprintf('\nNon-normal priors specified. %s with penalty uses a Laplace type of approximation.\n',options_mom_.mom.mom_method)
        fprintf('Only the prior mean and standard deviation are relevant, all other shape information, except for the parameter bounds, is ignored.\n\n')
        non_normal_priors=bayestopt_laplace.pshape~=3;
        bayestopt_laplace.pshape(non_normal_priors) = 3;
        bayestopt_laplace.p3(non_normal_priors) = -Inf*ones(sum(non_normal_priors),1);
        bayestopt_laplace.p4(non_normal_priors) = Inf*ones(sum(non_normal_priors),1);
        bayestopt_laplace.p6(non_normal_priors) = bayestopt_laplace.p1(non_normal_priors);
        bayestopt_laplace.p7(non_normal_priors) = bayestopt_laplace.p2(non_normal_priors);
        bayestopt_laplace.p5(non_normal_priors) = bayestopt_laplace.p1(non_normal_priors);
    end
    if any(isinf(bayestopt_laplace.p2)) %find infinite variance priors
        inf_var_pars=bayestopt_laplace.name(isinf(bayestopt_laplace.p2));
        disp_string=[inf_var_pars{1,:}];
        for ii=2:size(inf_var_pars,1)
            disp_string=[disp_string,', ',inf_var_pars{ii,:}];
        end
        fprintf('The parameter(s) %s have infinite prior variance. This implies a flat prior\n',disp_string)
        fprintf('Dynare disables the matrix singularity warning\n')
        if isoctave
            warning('off','Octave:singular-matrix');
        else
            warning('off','MATLAB:singularMatrix');
        end
    end
end

% Check for calibrated covariances before updating parameters
estim_params_ = check_for_calibrated_covariances(xparam0,estim_params_,M_);

% Checks on parameter calibration and initialization
xparam1_calib = get_all_parameters(estim_params_,M_); %get calibrated parameters
if ~any(isnan(xparam1_calib)) %all estimated parameters are calibrated
    estim_params_.full_calibration_detected=1;
else
    estim_params_.full_calibration_detected=0;
end
if options_mom_.use_calibration_initialization %set calibration as starting values
    if ~isempty(bayestopt_laplace) && all(bayestopt_laplace.pshape==0) && any(all(isnan([xparam1_calib xparam0]),2))
        error('method_of_moments: When using the use_calibration option with %s without prior, the parameters must be explicitly initialized.',options_mom_.mom.mom_method)
    else
        [xparam0,estim_params_]=do_parameter_initialization(estim_params_,xparam1_calib,xparam0); %get explicitly initialized parameters that have precedence over calibrated values
    end
end

% Check initialization
if ~isempty(bayestopt_laplace) && all(bayestopt_laplace.pshape==0) && any(isnan(xparam0))
    error('method_of_moments: %s without penalty requires all estimated parameters to be initialized, either in an estimated_params or estimated_params_init-block ',options_mom_.mom.mom_method)
end

% Set and check parameter bounds
if ~isempty(bayestopt_laplace) && any(bayestopt_laplace.pshape > 0)
    % Plot prior densities
    if ~options_mom_.nograph && options_mom_.plot_priors
        plot_priors(bayestopt_,M_,estim_params_,options_mom_)
        plot_priors(bayestopt_laplace,M_,estim_params_,options_mom_,'Laplace approximated priors')
    end
    % Set prior bounds
    Bounds = prior_bounds(bayestopt_laplace, options_mom_.prior_trunc);
    Bounds.lb = max(Bounds.lb,lb);
    Bounds.ub = min(Bounds.ub,ub);
else  % estimated parameters but no declared priors
    % No priors are declared so Dynare will estimate the parameters
    % with inequality constraints for the parameters.
    Bounds.lb = lb;
    Bounds.ub = ub;
    if options_mom_.mom.penalized_estimator
        fprintf('Penalized estimation turned off as you did not declare priors\n')
        options_mom_.mom.penalized_estimator = 0;
    end    
end
% Set correct bounds for standard deviations and corrrelations
param_of_interest=(1:length(xparam0))'<=estim_params_.nvx+estim_params_.nvn;
LB_below_0=(Bounds.lb<0 & param_of_interest);
Bounds.lb(LB_below_0)=0;
param_of_interest=(1:length(xparam0))'> estim_params_.nvx+estim_params_.nvn & (1:length(xparam0))'<estim_params_.nvx+estim_params_.nvn +estim_params_.ncx + estim_params_.ncn;
LB_below_minus_1=(Bounds.lb<-1 & param_of_interest);
UB_above_1=(Bounds.ub>1 & param_of_interest);
Bounds.lb(LB_below_minus_1)=-1; 
Bounds.ub(UB_above_1)=1; 

clear('bayestopt_','LB_below_0','LB_below_minus_1','UB_above_1','param_of_interest');%make sure stale structure cannot be used

% Test if initial values of the estimated parameters are all between the prior lower and upper bounds
if options_mom_.use_calibration_initialization
    try
        check_prior_bounds(xparam0,Bounds,M_,estim_params_,options_mom_,bayestopt_laplace)
    catch last_error
        fprintf('Cannot use parameter values from calibration as they violate the prior bounds.')
        rethrow(last_error);
    end
else
    check_prior_bounds(xparam0,Bounds,M_,estim_params_,options_mom_,bayestopt_laplace)
end

estim_params_= get_matrix_entries_for_psd_check(M_,estim_params_);

% Set sigma_e_is_diagonal flag (needed if the shocks block is not declared in the mod file).
M_.sigma_e_is_diagonal = true;
if estim_params_.ncx || any(nnz(tril(M_.Correlation_matrix,-1))) || isfield(estim_params_,'calibrated_covariances')
    M_.sigma_e_is_diagonal = false;
end

% storing prior parameters in MoM info structure for penalized minimization
oo_.prior.mean = bayestopt_laplace.p1;
oo_.prior.variance = diag(bayestopt_laplace.p2.^2);

% Set all parameters
M_ = set_all_parameters(xparam0,estim_params_,M_);

%provide warning if there is NaN in parameters
test_for_deep_parameters_calibration(M_);

% -------------------------------------------------------------------------
% Step 4: Checks and transformations for data
% -------------------------------------------------------------------------

% Check if datafile has same name as mod file
[~,name,~] = fileparts(options_mom_.datafile);
if strcmp(name,M_.fname)
    error('method_of_moments: Data-file and mod-file are not allowed to have the same name. Please change the name of the data file.')
end

% Build dataset
dataset_ = makedataset(options_mom_);

% set options for old interface from the ones for new interface
if ~isempty(dataset_)
    options_mom_.nobs = dataset_.nobs;
end

% Check length of data for estimation of second moments
if options_mom_.ar > options_mom_.nobs+1
    error('method_of_moments: Data set is too short to compute second moments');
end

% Provide info on data moments handling
fprintf('Computing data moments. Note that NaN values in the moments (due to leads and lags or missing data) are replaced by the mean of the corresponding moment\n');

% Get data moments for the method of moments
[oo_.mom.data_moments, oo_.mom.m_data] = mom.data_moments(dataset_.data, oo_, M_.matched_moments, options_mom_);

% Get shock series for SMM and set variance correction factor
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
[oo_.steady_state, new_steady_params, info] = evaluate_steady_state(oo_.steady_state,M_,options_mom_,oo_,steadystate_check_flag);
if info(1)
    fprintf('\nmethod_of_moments: The steady state at the initial parameters cannot be computed.\n')
    print_info(info, 0, options_mom_);
end

% check whether steady state file changes estimated parameters
if isfield(estim_params_,'param_vals') && ~isempty(estim_params_.param_vals)
    Model_par_varied=M_; %store M_ structure
    
    Model_par_varied.params(estim_params_.param_vals(:,1))=Model_par_varied.params(estim_params_.param_vals(:,1))*1.01; %vary parameters
    [~, new_steady_params_2] = evaluate_steady_state(oo_.steady_state,Model_par_varied,options_mom_,oo_,1);
    
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
if all(abs(oo_.steady_state(oo_.dr.order_var(oo_.dr.obs_var)))<1e-9)
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

for i = 1:length(optimizer_vec)
    if i == 1
        str = '- optimizer (mode_compute';
    else
        str = '            (additional_optimizer_steps';
    end
    switch optimizer_vec{i}
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
            if ischar(optimizer_vec{i})
                fprintf('\n  %s=%s): user-defined',str,optimizer_vec{i});
            else
                error('method_of_moments: Unknown optimizer, please contact the developers ')
            end
    end
    if options_mom_.silent_optimizer
        fprintf(' (silent)');
    end
    if strcmp(options_mom_.mom.mom_method,'GMM') && options_mom_.mom.analytic_jacobian && ismember(optimizer_vec{i},analytic_jacobian_optimizers)
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

    for optim_iter= 1:length(optimizer_vec)
        options_mom_.current_optimizer = optimizer_vec{optim_iter};
        if optimizer_vec{optim_iter}==0
            xparam1=xparam0; %no minimization, evaluate objective at current values
            fval = feval(objective_function, xparam1, Bounds, oo_, estim_params_, M_, options_mom_);
        else
            if optimizer_vec{optim_iter}==13
                options_mom_.vector_output = true;                
            else
                options_mom_.vector_output = false;                
            end
            if strcmp(options_mom_.mom.mom_method,'GMM') && options_mom_.mom.analytic_jacobian && ismember(optimizer_vec{optim_iter},analytic_jacobian_optimizers) %do this only for gradient-based optimizers
                options_mom_.mom.compute_derivs = true;
            else
                options_mom_.mom.compute_derivs = false;
            end
            
            [xparam1, fval, exitflag] = dynare_minimize_objective(objective_function, xparam0, optimizer_vec{optim_iter}, options_mom_, [Bounds.lb Bounds.ub], bayestopt_laplace.name, bayestopt_laplace, [],...
                                                                  Bounds, oo_, estim_params_, M_, options_mom_);
            if options_mom_.vector_output
                fval = fval'*fval;
            end
        end
        fprintf('\nStage %d Iteration %d: value of minimized moment distance objective function: %12.10f.\n',stage_iter,optim_iter,fval)
        if options_mom_.mom.verbose
            oo_.mom=display_estimation_results_table(xparam1,NaN(size(xparam1)),M_,options_mom_,estim_params_,bayestopt_laplace,oo_.mom,prior_dist_names,sprintf('%s (STAGE %d ITERATION %d) VERBOSE',options_mom_.mom.mom_method,stage_iter,optim_iter),sprintf('verbose_%s_stage_%d_iter_%d',lower(options_mom_.mom.mom_method),stage_iter,optim_iter));
        end
        xparam0=xparam1;
    end
    options_mom_.vector_output = false;    
    % Update M_ and DynareResults (in particular to get oo_.mom.model_moments)    
    M_ = set_all_parameters(xparam1,estim_params_,M_);
    if strcmp(options_mom_.mom.mom_method,'GMM') && options_mom_.mom.analytic_standard_errors
        options_mom_.mom.compute_derivs = true; % for GMM we compute derivatives analytically in the objective function with this flag        
    end
    [fval, ~, ~,~,~, oo_] = feval(objective_function, xparam1, Bounds, oo_, estim_params_, M_, options_mom_);
    options_mom_.mom.compute_derivs = false; % reset to not compute derivatives in objective function during optimization
    
    SE = mom.standard_errors(xparam1, objective_function, Bounds, oo_, estim_params_, M_, options_mom_, Woptflag);
    
    % Store results in output structure
    oo_.mom = display_estimation_results_table(xparam1,SE,M_,options_mom_,estim_params_,bayestopt_laplace,oo_.mom,prior_dist_names,sprintf('%s (STAGE %u)',options_mom_.mom.mom_method,stage_iter),sprintf('%s_stage_%u',lower(options_mom_.mom.mom_method),stage_iter));
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
    mom.check_plot(objective_function,xparam1,SE,options_mom_,M_,estim_params_,Bounds,bayestopt_laplace,...
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
