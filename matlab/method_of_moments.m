function [DynareResults, OptionsMoM] = method_of_moments(BayesInfo, DynareOptions, DynareResults, EstimatedParameters, Model, MatchedMoments, OptionsMoM)
%function [oo_, options_mom] = method_of_moments(M, options, oo, bayestopt, estim_params, matched_moments, options_mom)
% -------------------------------------------------------------------------
% This function performs a method of moments estimation with the following steps:
%   Step 0: Check if required structures and options exist
%   Step 1: - Prepare OptionsMoM structure
%           - Carry over Options from the preprocessor
%           - Other options that need to be initialized
%           - Get variable orderings and state space representation
%   Step 2: Checks and transformations for matched moments structure (preliminary)
%   Step 3: Checks and transformations for estimated parameters, priors, and bounds
%   Step 4: Checks and transformations for data
%   Step 5: checks for steady state at initial parameters
%   Step 6: checks for objective function at initial parameters
%   Step 7: Method of moments estimation: print some info, first-stage, and second-stage
%   Step 8: Clean up
% -------------------------------------------------------------------------
% This function is inspired by replication codes accompanied to the following papers:
%  o Andreasen, Fernández-Villaverde, Rubio-Ramírez (2018): "The Pruned State-Space System for Non-Linear DSGE Models: Theory and Empirical Applications", Review of Economic Studies, 85(1):1-49.
%  o Born, Pfeifer (2014): "Risk Matters: Comment", American Economic Review, 104(12):4231-4239.
%  o Mutschler (2018): "Higher-order statistics for DSGE models", Econometrics and Statistics, 6:44-56.
% =========================================================================
% INPUTS
%  o BayesInfo:              [structure] information about priors (bayestopt_)
%  o DynareOptions:          [structure] information about global options (options_)
%  o DynareResults:          [structure] storage for results (oo_)
%  o EstimatedParameters:    [structure] information about estimated parameters (estim_params_)
%  o Model:                  [structure] information about model (M_)
%  o MatchedMoments:         [structure] information about selected moments to match in estimation (matched_moments_)
%  o OptionsMoM:             [structure] information about settings specified by the user (options_mom_)
% -------------------------------------------------------------------------
% OUTPUTS
%  o DynareResults:          [structure] storage for results (oo_)
%  o OptionsMoM:             [structure] information about all used settings used in estimation (options_mom_)
% -------------------------------------------------------------------------
% This function is called by
%  o driver.m
% -------------------------------------------------------------------------
% This function calls
%  o check_for_calibrated_covariances.m
%  o check_prior_bounds.m
%  o do_parameter_initialization.m
%  o dynare_minimize.m
%  o evaluate_steady_state
%  o get_all_parameters.m
%  o makedataset.m
%  o method_of_moments_datamoments.m
%  o plot_priors.m
%  o print_info.m
%  o prior_bounds.m
%  o set_default_option.m
%  o set_prior.m
%  o set_state_space.m
%  o set_all_parameters
%  o list_of_parameters_calibrated_as_Inf
%  o list_of_parameters_calibrated_as_NaN
% =========================================================================
% Copyright (C) 2020 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------
% Author(s): 
% o Willi Mutschler (willi@mutschler.eu)
% o Johannes Pfeifer (jpfeifer@uni-koeln.de)
% =========================================================================

%% TO DO LIST
% - [ ] penalized estimation: how does penalized_estimator work? What is
%       penalized_estimator? Not all optimizer make use of this...what is special
%       about mode_compute=1 in objective functions. Do we need global objective_function_penalty_base in objective function
% - [ ] make csminwel work
% - [ ] why does lsqnonlin take less time in Andreasen toolbox?
% - [ ] recheck different optimizers if they are useful
% - [ ] test prefilter option (i.e. centered moments)
% - [ ] how to deal with logged_steady_state
% - [ ] mode_check plots?
% - [ ] test prior restrictions file
% - [ ] test user-specified weightning matrix
% - [ ] test non-symetric Sigma_e
% - [ ] which qz_criterium value?
% - [ ] are missing observations a problem? in method_of_moments_datamoments.m nan are replaced by mean of moment
% - [ ] check negative priors on std errors and above 1 for correlations
% - [ ] add measurement errors
% - [ ] add IRF matching
% - [ ] what are trend_coeffs and how to deal with them?
% - [ ] once interface is established: provide code to remove duplicate declared moments in matched moments block
% - [ ] test estimated_params_init(use calibration)
% - [ ] test estimated_params_bounds block
% - [ ] test what happens if all parameters will be estimated but some/all are not calibrated
% - [ ] speed up lyapunov equation by using doubling with old initial values
% - [ ] check what happens if parameters are set to INF or NAN
% - [ ] provide a table with dataMoments and final modelMoments 
% - [ ] check smm at order > 3 without pruning
% - [ ] provide option to use analytical derivatives to compute std errors (similar to what we already do in identification)
% - [ ] add Bayesian GMM/SMM estimation
% - [ ] useautocorr
% - [ ] instead of mom_steps, iterate over optimizers using additional_optimizer_steps

% -------------------------------------------------------------------------
% Step 0: Check if required structures and options exist
% -------------------------------------------------------------------------
if isempty(EstimatedParameters) % structure storing the info about estimated parameters in the estimated_params block
    error('method_of_moments: You need to provide an ''estimated_params'' block')
end
if isempty(MatchedMoments) % structure storing the moments used for the method of moments estimation
    error('method_of_moments: You need to provide a ''matched_moments'' block')
end
if ~isempty(BayesInfo) && any(BayesInfo.pshape==0) && any(BayesInfo.pshape~=0)
    error('method_of_moments: Estimation must be either fully classical or fully Bayesian. Maybe you forgot to specify a prior distribution.')
end
fprintf('\n==== Method of Moments Estimation ====\n\n')

% -------------------------------------------------------------------------
% Step 1a: Prepare OptionsMoM structure
% -------------------------------------------------------------------------
% OptionsMoM is local and contains default and user-specified values for 
% all settings needed for the method of moments estimation. Some options,
% though, are set by the preprocessor into options_ and we copy these over.
% The idea is to be independent of options_ and have full control of the
% estimation instead of possibly having to deal with options chosen somewhere
% else in the mod file.

% Method of Moments estimation options that can be set by the user in the mod file, otherwise default values are provided
if strcmp(OptionsMoM.mom_method,'GMM') || strcmp(OptionsMoM.mom_method,'SMM')
    OptionsMoM = set_default_option(OptionsMoM,'bartlett_kernel_lag',20);             % bandwith in optimal weighting matrix
    OptionsMoM = set_default_option(OptionsMoM,'order',1);                            % order of Taylor approximation in perturbation
    OptionsMoM = set_default_option(OptionsMoM,'penalized_estimator',false);          % @wmutschl: provide description
    OptionsMoM = set_default_option(OptionsMoM,'pruning',true);                       % use pruned state space system at higher-order
    OptionsMoM = set_default_option(OptionsMoM,'verbose',false);                       % display and store intermediate estimation results
    OptionsMoM = set_default_option(OptionsMoM,'weighting_matrix','identity_matrix'); % weighting matrix in moments distance objective function; possible values: OPTIMAL|IDENTITY_MATRIX|DIAGONAL|filename
    OptionsMoM = set_default_option(OptionsMoM,'additional_optimizer_steps',[]);                    % Number of iterations in the steps of the 2-step feasible method of moments
    % Checks for perturbation order
    if OptionsMoM.order < 1
        error('method_of_moments:: The order of the Taylor approximation cannot be 0!')
    elseif OptionsMoM.order > 2 && (~isfield(OptionsMoM,'k_order_solver') || ~OptionsMoM.k_order_solver)
        fprintf('method_of_moments: For perturbation order k>2, we add the k_order_solver option.\n');
        OptionsMoM.k_order_solver = 1;
    end
end

if strcmp(OptionsMoM.mom_method,'SMM')
    objective_function = str2func('method_of_moments_SMM');
    OptionsMoM = set_default_option(OptionsMoM,'bounded_shock_support',false);        % trim shocks in simulation to +- 2 stdev
    OptionsMoM = set_default_option(OptionsMoM,'drop',500);                           % number of periods dropped at beginning of simulation
    OptionsMoM = set_default_option(OptionsMoM,'seed',24051986);                      % seed used in simulations
    OptionsMoM = set_default_option(OptionsMoM,'simulation_multiple',5);              % multiple of the data length used for simulation
    if OptionsMoM.simulation_multiple < 1
        fprintf('The simulation horizon is shorter than the data. Dynare resets the simulation_multiple to 2.\n')
        OptionsMoM.smm.simulation_multiple = 2;
    end
end
if strcmp(OptionsMoM.mom_method,'GMM')
    objective_function = str2func('method_of_moments_GMM');
    % Check for pruning with GMM at higher order
    if OptionsMoM.order > 1 && ~OptionsMoM.pruning
        fprintf('GMM at higher order only works with pruning, so we set pruning option to 1.\n');
        OptionsMoM.pruning = true;
    end
end

% General options that can be set by the user in the mod file, otherwise default values are provided
OptionsMoM = set_default_option(OptionsMoM,'dirname',Model.fname); % directory in which to store estimation output
OptionsMoM = set_default_option(OptionsMoM,'graph_format','eps');  % specify the file format(s) for graphs saved to disk
OptionsMoM = set_default_option(OptionsMoM,'nodisplay',false);     % do not display the graphs, but still save them to disk
OptionsMoM = set_default_option(OptionsMoM,'nograph',false);       % do not create graphs (which implies that they are not saved to the disk nor displayed)
OptionsMoM = set_default_option(OptionsMoM,'noprint',false);       % do not print output to console
OptionsMoM = set_default_option(OptionsMoM,'plot_priors',true);    % control plotting of priors
OptionsMoM = set_default_option(OptionsMoM,'prior_trunc',1e-10);   % probability of extreme values of the prior density that is ignored when computing bounds for the parameters
OptionsMoM = set_default_option(OptionsMoM,'TeX',false);           % print TeX tables and graphics

% Data and model options that can be set by the user in the mod file, otherwise default values are provided
OptionsMoM = set_default_option(OptionsMoM,'first_obs',1);     % number of first observation
OptionsMoM = set_default_option(OptionsMoM,'logdata',false);   % if loglinear is set, this option is necessary if the user provides data already in logs, otherwise the log transformation will be applied twice (this may result in complex data)
OptionsMoM = set_default_option(OptionsMoM,'loglinear',false); % computes a log-linear approximation of the model instead of a linear approximation
OptionsMoM = set_default_option(OptionsMoM,'nobs',NaN);        % number of observations
OptionsMoM = set_default_option(OptionsMoM,'prefilter',false); % demean each data series by its empirical mean and use centered moments
OptionsMoM = set_default_option(OptionsMoM,'xls_sheet',1);     % name of sheet with data in Excel
OptionsMoM = set_default_option(OptionsMoM,'xls_range','');    % range of data in Excel sheet
% Recursive estimation and forecast are not supported
if numel(OptionsMoM.nobs)>1
    error('method_of_moments: Recursive estimation and forecast for samples is not supported. Please set an integer as ''nobs''.');
end
if numel(OptionsMoM.first_obs)>1
    error('method_of_moments: Recursive estimation and forecast for samples is not supported. Please set an integer as ''first_obs''.');
end

% Optimization options that can be set by the user in the mod file, otherwise default values are provided
if strcmp(OptionsMoM.mom_method, 'GMM')
    OptionsMoM = set_default_option(OptionsMoM,'analytic_derivation',0); % use analytic derivatives to compute standard errors for GMM
elseif isfield(OptionsMoM,'analytic_derivation')
    fprintf('Only GMM supports analytic derivation to compute standard errors, we reset ''analytic_derivation'' to 0.\n')
    OptionsMoM.analytic_derivation = 0;
else
    OptionsMoM.analytic_derivation = 0;
end
OptionsMoM = set_default_option(OptionsMoM,'huge_number',1e7);        % value for replacing the infinite bounds on parameters by finite numbers. Used by some optimizers for numerical reasons
OptionsMoM = set_default_option(OptionsMoM,'mode_compute',13);        % specifies the optimizer for minimization of moments distance
OptionsMoM = set_default_option(OptionsMoM,'optim_opt',[]);           % a list of NAME and VALUE pairs to set options for the optimization routines. Available options depend on mode_compute
OptionsMoM = set_default_option(OptionsMoM,'silent_optimizer',false); % run minimization of moments distance silently without displaying results or saving files in between
if ~isfield(OptionsMoM,'dynatol')
    OptionsMoM.dynatol = {};
end
OptionsMoM.dynatol = set_default_option(OptionsMoM.dynatol,'f', 1e-5);% convergence criterion on function value for numerical differentiation
OptionsMoM.dynatol = set_default_option(OptionsMoM.dynatol,'x', 1e-5);% convergence criterion on funciton input for numerical differentiation

% Numerical algorithms options that can be set by the user in the mod file, otherwise default values are provided
OptionsMoM = set_default_option(OptionsMoM,'aim_solver',false);                     % Use AIM algorithm to compute perturbation approximation
OptionsMoM = set_default_option(OptionsMoM,'dr_cycle_reduction',false);             % use cycle reduction algorithm to solve the polynomial equation for retrieving the coefficients associated to the endogenous variables in the decision rule
OptionsMoM = set_default_option(OptionsMoM,'dr_cycle_reduction_tol',1e-7);          % convergence criterion used in the cycle reduction algorithm
OptionsMoM = set_default_option(OptionsMoM,'dr_logarithmic_reduction',false);       % use logarithmic reduction algorithm to solve the polynomial equation for retrieving the coefficients associated to the endogenous variables in the decision rule
OptionsMoM = set_default_option(OptionsMoM,'dr_logarithmic_reduction_maxiter',100); % maximum number of iterations used in the logarithmic reduction algorithm
OptionsMoM = set_default_option(OptionsMoM,'dr_logarithmic_reduction_tol',1e-12);   % convergence criterion used in the cycle reduction algorithm
OptionsMoM = set_default_option(OptionsMoM,'lyapunov_db',false);                    % doubling algorithm (disclyap_fast) to solve Lyapunov equation to compute variance-covariance matrix of state variables
OptionsMoM = set_default_option(OptionsMoM,'lyapunov_fp',false);                    % fixed-point algorithm to solve Lyapunov equation to compute variance-covariance matrix of state variables
OptionsMoM = set_default_option(OptionsMoM,'lyapunov_srs',false);                   % square-root-solver (dlyapchol) algorithm to solve Lyapunov equation to compute variance-covariance matrix of state variables
OptionsMoM = set_default_option(OptionsMoM,'lyapunov_complex_threshold',1e-15);     % complex block threshold for the upper triangular matrix in symmetric Lyapunov equation solver
OptionsMoM = set_default_option(OptionsMoM,'lyapunov_fixed_point_tol',1e-10);       % convergence criterion used in the fixed point Lyapunov solver
OptionsMoM = set_default_option(OptionsMoM,'lyapunov_doubling_tol',1e-16);          % convergence criterion used in the doubling algorithm
OptionsMoM = set_default_option(OptionsMoM,'sylvester_fp',false);                   % determines whether to use fixed point algorihtm to solve Sylvester equation (gensylv_fp), faster for large scale models
OptionsMoM = set_default_option(OptionsMoM,'sylvester_fixed_point_tol',1e-12);      % convergence criterion used in the fixed point Sylvester solver
OptionsMoM = set_default_option(OptionsMoM,'qz_criterium',1-1e-6);                  % value used to split stable from unstable eigenvalues in reordering the Generalized Schur decomposition used for solving first order problems [IS THIS CORRET @wmutschl]
OptionsMoM = set_default_option(OptionsMoM,'qz_zero_threshold',1e-6);               % value used to test if a generalized eigenvalue is 0/0 in the generalized Schur decomposition

% -------------------------------------------------------------------------
% Step 1b: Options that are set by the preprocessor and (probably) need to be carried over
% -------------------------------------------------------------------------

% options related to VAROBS
if ~isfield(DynareOptions,'varobs')
    error('method_of_moments: VAROBS statement is missing!')
else
    OptionsMoM.varobs  = DynareOptions.varobs;     % observable variables in declaration order
    OptionsMoM.obs_nbr = length(OptionsMoM.varobs);% number of observed variables
    % Check that each declared observed variable is also an endogenous variable
    for i = 1:OptionsMoM.obs_nbr
        if ~any(strcmp(OptionsMoM.varobs{i},Model.endo_names))
            error(['method_of_moments: Unknown variable (' OptionsMoM.varobs{i} ')!'])
        end
    end

    % Check that a variable is not declared as observed more than once.
    if length(unique(OptionsMoM.varobs))<length(OptionsMoM.varobs)
        for i = 1:OptionsMoM.obs_nbr
            if sum(strcmp(OptionsMoM.varobs{i},OptionsMoM.varobs))>1
                error(['method_of_moments: A variable cannot be declared as observed more than once (' OptionsMoM.varobs{i} ')!'])
            end
        end
    end
end

% options related to variable declarations
if ~isfield(DynareOptions,'trend_coeffs')
    %BayesInfo.with_trend = 0;
else
    error('method_of_moments: %s does not allow for trend in data',OptionsMoM.mom_method)
end

% options related to estimated_params and estimated_params_init
OptionsMoM.use_calibration_initialization = DynareOptions.use_calibration_initialization;

% options related to model; block
OptionsMoM.linear   = DynareOptions.linear;
OptionsMoM.use_dll  = DynareOptions.use_dll;
OptionsMoM.block    = DynareOptions.block;
OptionsMoM.bytecode = DynareOptions.bytecode;

% options related to steady; command
OptionsMoM.homotopy_force_continue = DynareOptions.homotopy_force_continue;
OptionsMoM.homotopy_mode           = DynareOptions.homotopy_mode;
OptionsMoM.homotopy_steps          = DynareOptions.homotopy_steps;
OptionsMoM.logged_steady_state     = DynareOptions.logged_steady_state; % @wmutschl: when and how does this get set?
OptionsMoM.markowitz               = DynareOptions.markowitz;
OptionsMoM.solve_algo              = DynareOptions.solve_algo;
OptionsMoM.solve_tolf              = DynareOptions.solve_tolf;
OptionsMoM.steady                  = DynareOptions.steady;
OptionsMoM.steadystate             = DynareOptions.steadystate;
OptionsMoM.steadystate_flag        = DynareOptions.steadystate_flag;

% options related to dataset
OptionsMoM.dataset        = DynareOptions.dataset;
OptionsMoM.initial_period = DynareOptions.initial_period;

% options related to endogenous prior restrictions
OptionsMoM.endogenous_prior_restrictions.irf    = {};
OptionsMoM.endogenous_prior_restrictions.moment = {};
if ~isempty(DynareOptions.endogenous_prior_restrictions.irf) && ~isempty(DynareOptions.endogenous_prior_restrictions.moment)
    fprintf('Endogenous prior restrictions are not supported yet and will be skipped.\n')
end


% -------------------------------------------------------------------------
% Step 1c: Options related to optimizers
% -------------------------------------------------------------------------
% mode_compute = 2
OptionsMoM.saopt            = DynareOptions.saopt;
% mode_compute = 4
OptionsMoM.csminwel         = DynareOptions.csminwel;
if OptionsMoM.mode_compute == 4
    error('method_of_moments:optimizer','method_of_moments: csminwel optimizer (mode_compute=4) is not yet supported (due to penalized_objective handling).\n                   Choose a different optimizer, e.g. lsqnonlin (mode_compute=13), fminsearch (mode_compute=7), SolveOpt (mode_compute=101).')
end
% mode_compute = 5 (not yet)
if OptionsMoM.mode_compute == 5
    error('method_of_moments:optimizer','method_of_moments: newrat optimizer (mode_compute=5) is not yet supported.\n                   Choose a different optimizer, e.g. lsqnonlin (mode_compute=13), fminsearch (mode_compute=7), SolveOpt (mode_compute=101).')
end
% mode_compute = 6
if OptionsMoM.mode_compute == 6
    error('method_of_moments:optimizer','method_of_moments: mode_compute=6 is not compatible with a method of moments estimation.\n                   Choose a different optimizer, e.g. lsqnonlin (mode_compute=13), fminsearch (mode_compute=7), SolveOpt (mode_compute=101).')
end
% mode_compute = 8
OptionsMoM.simplex          = DynareOptions.simplex;
% mode_compute = 9
OptionsMoM.cmaes            = DynareOptions.cmaes;
% mode_compute = 10
OptionsMoM.simpsa           = DynareOptions.simpsa;
% mode_compute = 11
if OptionsMoM.mode_compute == 11
    error('method_of_moments:optimizer','method_of_moments: mode_compute=11 is not compatible with a method of moments estimation.\n                   Choose a different optimizer, e.g. lsqnonlin (mode_compute=13), fminsearch (mode_compute=7), SolveOpt (mode_compute=101).')
end
% mode_compute = 12
OptionsMoM.particleswarm    = DynareOptions.particleswarm;
if OptionsMoM.mode_compute == 12
    error('method_of_moments:optimizer','method_of_moments: mode_compute=12 is not yet supported (due to penalized_objective handling).\n                   Choose a different optimizer, e.g. lsqnonlin (mode_compute=13), fminsearch (mode_compute=7), SolveOpt (mode_compute=101).')
end
% mode_compute = 101
OptionsMoM.solveopt         = DynareOptions.solveopt;

OptionsMoM.gradient_method  = DynareOptions.gradient_method;
OptionsMoM.gradient_epsilon = DynareOptions.gradient_epsilon;

% -------------------------------------------------------------------------
% Step 1d: Other options that need to be initialized
% -------------------------------------------------------------------------
OptionsMoM.initialize_estimated_parameters_with_the_prior_mode = 0; % needed by set_prior.m
OptionsMoM.figures.textwidth = 0.8; %needed by plot_priors.m
OptionsMoM.ramsey_policy = 0; % needed by evaluate_steady_state
OptionsMoM.debug = false; %needed by resol.m
OptionsMoM.risky_steadystate = false; %needed by resol
OptionsMoM.threads = DynareOptions.threads; %needed by resol
OptionsMoM.jacobian_flag = true;
OptionsMoM.gstep = DynareOptions.gstep;
OptionsMoM.solve_tolf = DynareOptions.solve_tolf;
OptionsMoM.solve_tolx = DynareOptions.solve_tolx;


% options_mom.dsge_var          = false; %needed by check_list_of_variables
% options_mom.bayesian_irf      = false; %needed by check_list_of_variables
% options_mom.moments_varendo   = false; %needed by check_list_of_variables
% options_mom.smoother          = false; %needed by check_list_of_variables
% options_mom.filter_step_ahead = [];  %needed by check_list_of_variables
% options_mom.forecast = 0;
%OptionsMoM = set_default_option(OptionsMoM,'endo_vars_for_moment_computations_in_estimation',[]);

% -------------------------------------------------------------------------
% Step 1e: Get variable orderings and state space representation
% -------------------------------------------------------------------------
DynareResults.dr = set_state_space(DynareResults.dr,Model,OptionsMoM);
% Get index of observed variables in DR order
DynareResults.dr.obs_var = [];
for i=1:OptionsMoM.obs_nbr
    DynareResults.dr.obs_var = [DynareResults.dr.obs_var; find(strcmp(OptionsMoM.varobs{i}, Model.endo_names(DynareResults.dr.order_var)))];
end

% -------------------------------------------------------------------------
% Step 2: Checks and transformations for matched moments structure (preliminary)
% -------------------------------------------------------------------------
% Note that we do not have a preprocessor interface yet for this, so this
% will need much improvement later on. @wmutschl
if strcmp(OptionsMoM.mom_method, 'GMM')
    % Initialize indices
    OptionsMoM.index.E_y       = false(OptionsMoM.obs_nbr,1);                    %unconditional first order product moments    
    OptionsMoM.index.E_yy      = false(OptionsMoM.obs_nbr,OptionsMoM.obs_nbr);   %unconditional second order product moments
    OptionsMoM.index.E_yyt     = false(OptionsMoM.obs_nbr,OptionsMoM.obs_nbr,0); %unconditional temporal second order product moments
    OptionsMoM.index.E_y_pos   = zeros(OptionsMoM.obs_nbr,1);                    %position in matched moments block
    OptionsMoM.index.E_yy_pos  = zeros(OptionsMoM.obs_nbr,OptionsMoM.obs_nbr);   %position in matched moments block
    OptionsMoM.index.E_yyt_pos = zeros(OptionsMoM.obs_nbr,OptionsMoM.obs_nbr,0); %position in matched moments block
end

for jm=1:size(MatchedMoments,1)
    % higher-order product moments not supported yet for GMM
    if strcmp(OptionsMoM.mom_method, 'GMM') && sum(MatchedMoments{jm,3}) > 2
        error('method_of_moments: GMM does not yet support product moments higher than 2. Change row %d in ''matched_moments'' block.',jm);
    end    
    % Check if declared variables are also observed (needed as otherwise the dataset variables won't coincide)
    if any(~ismember(DynareResults.dr.inv_order_var(MatchedMoments{jm,1})', DynareResults.dr.obs_var))
        error('method_of_moments: Variables in row %d in ''matched_moments'' block need to be declared as VAROBS.', jm)
    end
    
    if strcmp(OptionsMoM.mom_method, 'GMM')
    % Check (for now) that only lags are declared
        if any(MatchedMoments{jm,2}>0)
            error('method_of_moments: Leads in row %d in the ''matched_moments'' block are not supported for GMM, shift the moments and declare only lags.', jm)
        end
        % Check (for now) that first declared variable has zero lag
        if MatchedMoments{jm,2}(1)~=0
            error('method_of_moments: The first variable declared in row %d in the ''matched_moments'' block is not allowed to have a lead or lag for GMM;\n                   reorder the variables in the row such that the first variable has zero lag!',jm)
        end
        vars = DynareResults.dr.inv_order_var(MatchedMoments{jm,1})';        
        if sum(MatchedMoments{jm,3}) == 1
            % First-order product moment
            vpos = (DynareResults.dr.obs_var == vars);
            OptionsMoM.index.E_y(vpos,1) = true;
            OptionsMoM.index.E_y_pos(vpos,1) = jm;            
        elseif sum(MatchedMoments{jm,3}) == 2
            % Second-order product moment
            idx1 = (DynareResults.dr.obs_var == vars(1));
            idx2 = (DynareResults.dr.obs_var == vars(2));
            lag1 = MatchedMoments{jm,2}(1);
            lag2 = MatchedMoments{jm,2}(2);
            if lag1==0 && lag2==0 % contemporenous covariance matrix
                OptionsMoM.index.E_yy(idx1,idx2) = true;
                OptionsMoM.index.E_yy(idx2,idx1) = true;
                OptionsMoM.index.E_yy_pos(idx1,idx2) = jm;
                OptionsMoM.index.E_yy_pos(idx2,idx1) = jm;
            elseif lag1==0 && lag2 < 0
                OptionsMoM.index.E_yyt(idx1,idx2,-lag2) = true;
                OptionsMoM.index.E_yyt_pos(idx1,idx2,-lag2) = jm;
            end
        end
    end
end

% @wmutschl: add check for duplicate moments by using the cellfun and unique functions
if strcmp(OptionsMoM.mom_method,'GMM')
    %Remove duplicate elements
    UniqueMomIdx = [nonzeros(OptionsMoM.index.E_y_pos); nonzeros(triu(OptionsMoM.index.E_yy_pos)); nonzeros(OptionsMoM.index.E_yyt_pos)];
    DuplicateMoms = setdiff(1:size(MatchedMoments,1),UniqueMomIdx);
    if ~isempty(DuplicateMoms)
        fprintf('Found and removed duplicate declared moments in ''matched_moments'' block in rows: %s.\n',num2str(DuplicateMoms))
    end
    %reorder MatchedMoments to be compatible with OptionsMoM.index
    MatchedMoments = MatchedMoments(UniqueMomIdx,:);
else    
    fprintf('For SMM we do not check yet for duplicate moment declarations in ''matched_moments'' block. You need to check this manually.\n\n')
end

% Check if both prefilter and first moments were specified
OptionsMoM.first_moment_indicator = find(cellfun(@(x) sum(abs(x))==1,MatchedMoments(:,3)))';
if OptionsMoM.prefilter && ~isempty(OptionsMoM.first_moment_indicator)
    fprintf('Centered moments requested (prefilter option is set); therefore, ignore declared first moments in ''matched_moments'' block in rows: %s.\n',num2str(OptionsMoM.index.E_y_pos'));
    MatchedMoments = MatchedMoments(OptionsMoM.first_moment_indicator,:); %remove first moments entries
    OptionsMoM.first_moment_indicator = [];
end
OptionsMoM.mom_nbr = size(MatchedMoments,1);



% Get maximum lag number for autocovariances/autocorrelations
OptionsMoM.ar = max(cellfun(@max,MatchedMoments(:,2))) - min(cellfun(@min,MatchedMoments(:,2)));

% -------------------------------------------------------------------------
% Step 3: Checks and transformations for estimated parameters, priors, and bounds
% -------------------------------------------------------------------------

% Set priors over the estimated parameters
if ~isempty(EstimatedParameters) && ~(isfield(EstimatedParameters,'nvx') && (size(EstimatedParameters.var_exo,1)+size(EstimatedParameters.var_endo,1)+size(EstimatedParameters.corrx,1)+size(EstimatedParameters.corrn,1)+size(EstimatedParameters.param_vals,1))==0)
    [xparam0, EstimatedParameters, BayesInfo, lb, ub, Model] = set_prior(EstimatedParameters, Model, OptionsMoM);
end

% Check measurement errors
if EstimatedParameters.nvn || EstimatedParameters.ncn
    error('method_of_moments: moment estimation does not support measurement error(s) yet. Please specifiy them as a structural shock.')
end

% Check if a _prior_restrictions.m file exists
if exist([Model.fname '_prior_restrictions.m'])
    OptionsMoM.prior_restrictions.status = 1;
    OptionsMoM.prior_restrictions.routine = str2func([Model.fname '_prior_restrictions']);
end

% Check on specified priors and penalized estimation
if ~isempty(EstimatedParameters) && ~(isfield(EstimatedParameters,'nvx') && (size(EstimatedParameters.var_exo,1)+size(EstimatedParameters.var_endo,1)+size(EstimatedParameters.corrx,1)+size(EstimatedParameters.corrn,1)+size(EstimatedParameters.param_vals,1))==0)
    if any(BayesInfo.pshape > 0) % prior specified
        if any(setdiff([0;BayesInfo.pshape],[0,3]))
            if ~OptionsMoM.penalized_estimator
                fprintf('\nPriors were specified, but the penalized_estimator-option was not set.\n')
                fprintf('Dynare sets penalized_estimator to 1. Conducting %s with penalty.\n',OptionsMoM.mom_method)
                OptionsMoM.penalized_estimator=1;
            end
            fprintf('\nNon-normal priors specified. %s with penalty uses a Laplace type of approximation.\n',OptionsMoM.mom_method)
            fprintf('Only the prior mean and standard deviation are relevant, all other shape information, except for the parameter bounds, is ignored.\n\n')
        end
        if any(isinf(BayesInfo.p2))
            inf_var_pars=BayesInfo.name(isinf(BayesInfo.p2));
            disp_string=[inf_var_pars{1,:}];
            for ii=2:size(inf_var_pars,1)
                disp_string=[disp_string,', ',inf_var_pars{ii,:}];
            end
            fprintf('The parameter(s) %s have infinite prior variance. This implies a flat prior\n',disp_string)
            fprintf('Dynare disables the matrix singularity warning\n')
            warning('off','MATLAB:singularMatrix');
        end
    end
end
if OptionsMoM.penalized_estimator && OptionsMoM.mode_compute==13
    error('method_of_moments: Penalized estimator option is not compatible with mode_compute=13. Deactivate penalized_estimator, don''t use priors, or choose a different optimizer.')
end


% Check for calibrated covariances before updating parameters
if ~isempty(EstimatedParameters) && ~(isfield(EstimatedParameters,'nvx') && sum(EstimatedParameters.nvx+EstimatedParameters.nvn+EstimatedParameters.ncx+EstimatedParameters.ncn+EstimatedParameters.np)==0)
    EstimatedParameters = check_for_calibrated_covariances(xparam0,EstimatedParameters,Model);
end

% Checks on parameter calibration and initialization
xparam1_calib = get_all_parameters(EstimatedParameters,Model); %get calibrated parameters
if ~any(isnan(xparam1_calib)) %all estimated parameters are calibrated
    EstimatedParameters.full_calibration_detected=1;
else
    EstimatedParameters.full_calibration_detected=0;
end
if OptionsMoM.use_calibration_initialization %set calibration as starting values
    if ~isempty(BayesInfo) && all(BayesInfo.pshape==0) && any(all(isnan([xparam1_calib xparam0]),2))
        error('method_of_moments: When using the use_calibration option with %s without prior, the parameters must be properly initialized.',OptionsMoM.mom_method)
    else
        [xparam0,EstimatedParameters]=do_parameter_initialization(EstimatedParameters,xparam1_calib,xparam0); %get explicitly initialized parameters that have precedence to calibrated values
    end
end

% Check on initialization @wmutschl: why without penalty?
if ~isempty(BayesInfo) && all(BayesInfo.pshape==0) && any(isnan(xparam0))
    error('method_of_moments: %s without penalty requires all estimated parameters to be initialized, either in an estimated_params or estimated_params_init-block ',OptionsMoM.mom_method)
end

% Set and check parameter bounds
if ~isempty(EstimatedParameters) && ~(all(strcmp(fieldnames(EstimatedParameters),'full_calibration_detected'))  || (isfield(EstimatedParameters,'nvx') && sum(EstimatedParameters.nvx+EstimatedParameters.nvn+EstimatedParameters.ncx+EstimatedParameters.ncn+EstimatedParameters.np)==0))
    if ~isempty(BayesInfo) && any(BayesInfo.pshape > 0)
        % Plot prior densities
        if ~OptionsMoM.nograph && OptionsMoM.plot_priors
            plot_priors(BayesInfo,Model,EstimatedParameters,OptionsMoM)
        end
        % Set prior bounds
        Bounds = prior_bounds(BayesInfo, OptionsMoM.prior_trunc);
        Bounds.lb = max(Bounds.lb,lb);
        Bounds.ub = min(Bounds.ub,ub);
    else  % estimated parameters but no declared priors
          % No priors are declared so Dynare will estimate the parameters
          % with inequality constraints for the parameters.
        Bounds.lb = lb;
        Bounds.ub = ub;
    end
    % Test if initial values of the estimated parameters are all between the prior lower and upper bounds
    if OptionsMoM.use_calibration_initialization
        try
            check_prior_bounds(xparam0,Bounds,Model,EstimatedParameters,OptionsMoM,BayesInfo)
        catch
            e = lasterror();
            fprintf('Cannot use parameter values from calibration as they violate the prior bounds.')
            rethrow(e);
        end
    else
        check_prior_bounds(xparam0,Bounds,Model,EstimatedParameters,OptionsMoM,BayesInfo)
    end
end

% Check symmetric Sigma_e
Sigma_e_non_zero_rows    = find(~all(Model.Sigma_e==0,1));
Sigma_e_non_zero_columns = find(~all(Model.Sigma_e==0,2));
if ~isequal(Sigma_e_non_zero_rows,Sigma_e_non_zero_columns')
    error('method_of_moments: Structual error matrix not symmetric')
end
if isfield(EstimatedParameters,'var_exo') && ~isempty(EstimatedParameters.var_exo)
    EstimatedParameters.Sigma_e_entries_to_check_for_positive_definiteness = union(Sigma_e_non_zero_rows,EstimatedParameters.var_exo(:,1));
else
    EstimatedParameters.Sigma_e_entries_to_check_for_positive_definiteness = Sigma_e_non_zero_rows;
end

% Set sigma_e_is_diagonal flag (needed if the shocks block is not declared in the mod file).
Model.sigma_e_is_diagonal = true;
if EstimatedParameters.ncx || any(nnz(tril(Model.Correlation_matrix,-1))) || isfield(EstimatedParameters,'calibrated_covariances')
    Model.sigma_e_is_diagonal = false;
end

% Provide some warnings on standard errors and correlations of shocks
if any(BayesInfo.pshape)        
    if EstimatedParameters.nvx && any(BayesInfo.p3(1:EstimatedParameters.nvx)<0)
        warning('Your prior allows for negative standard deviations for structural shocks. It is recommended to change your prior.')
    end
    offset=EstimatedParameters.nvx;        
    if EstimatedParameters.nvn && any(BayesInfo.p3(1+offset:offset+EstimatedParameters.nvn)<0)
        warning('Your prior allows for negative standard deviations for measurement error. It is recommended to change your prior.')
    end
    offset = EstimatedParameters.nvx+EstimatedParameters.nvn;        
    if EstimatedParameters.ncx && (any(BayesInfo.p3(1+offset:offset+EstimatedParameters.ncx)<-1) || any(BayesInfo.p4(1+offset:offset+EstimatedParameters.ncx)>1))
        warning('Your prior allows for correlations between structural shocks larger than +-1 and will not integrate to 1 due to truncation. Please change your prior')
    end
    offset = EstimatedParameters.nvx+EstimatedParameters.nvn+EstimatedParameters.ncx;        
    if EstimatedParameters.ncn && (any(BayesInfo.p3(1+offset:offset+EstimatedParameters.ncn)<-1) || any(BayesInfo.p4(1+offset:offset+EstimatedParameters.ncn)>1))
        warning('Your prior allows for correlations between measurement errors larger than +-1 and will not integrate to 1 due to truncation. Please change your prior')
    end
end   

% storing prior parameters in MoM info structure for penalized minimization
DynareResults.prior.pnames = {''; 'beta'; 'gamm'; 'norm'; 'invg'; 'unif'; 'invg2'; ''; 'weibl'};
DynareResults.prior.p1 = BayesInfo.p1;
DynareResults.prior.p2 = BayesInfo.p2;
% DynareResults.prior.mode                   = BayesInfo.p5;
% DynareResults.prior.variance               = diag(BayesInfo.p2.^2);
% DynareResults.prior.hyperparameters.first  = BayesInfo.p6;
% DynareResults.prior.hyperparameters.second = BayesInfo.p7;


% Set all parameters
Model = set_all_parameters(xparam0,EstimatedParameters,Model);


% -------------------------------------------------------------------------
% Step 4: Checks and transformations for data
% -------------------------------------------------------------------------

% Check if datafile has same name as mod file
if ~isempty(OptionsMoM.datafile)
    [~,name,~] = fileparts(OptionsMoM.datafile);
    if strcmp(name,Model.fname)
        error('method_of_moments: Data-file and mod-file are not allowed to have the same name. Please change the name of the data file.')
    end
end

% Build dataset
[DynareDataset, ~, ~] = makedataset(OptionsMoM);

% set options for old interface from the ones for new interface
if ~isempty(DynareDataset)
    OptionsMoM.nobs = DynareDataset.nobs;
end

% missing observations are not supported yet
if any(any(isnan(DynareDataset.data)))
    error('method_of_moments: missing observations are not supported')
end

% Check length of data for estimation of second moments
if OptionsMoM.ar > OptionsMoM.nobs+1
    error('method_of_moments: Data set is too short to compute second moments');
end

% Get data moments for the method of moments
[DynareResults.mom.dataMoments, DynareResults.mom.m_data] = method_of_moments_datamoments(DynareDataset.data, DynareResults, MatchedMoments, OptionsMoM);

if OptionsMoM.prefilter
    if sum(abs(DynareDataset.mean))/DynareDataset.nobs >1e-9
        fprintf('The mean of the data is:\n')
        disp(DynareDataset.mean);
        error('method_of_moments: You are trying to perform a method of moments estimation with centered moments (prefilter=1) using uncentered data.')
    end
elseif ~isempty(OptionsMoM.first_moment_indicator)    
    if sum(abs(DynareResults.mom.dataMoments(OptionsMoM.first_moment_indicator)))/sum(OptionsMoM.first_moment_indicator) <1e-2
        fprintf('The mean of the data for which Dynare tries to match first moments is:\n')
        disp(DynareResults.mom.dataMoments(OptionsMoM.first_moment_indicator)');
        warning('method_of_moments:data_mean_zero','method_of_moments: You are trying to perform a method of moments estimation with uncentered moments (prefilter=0),\n         but the data are (almost) mean 0. Check if this is desired.')
    end    
end

% Get shock series for SMM and set variance correction factor
if strcmp(OptionsMoM.mom_method,'SMM')
    OptionsMoM.long = round(OptionsMoM.simulation_multiple*OptionsMoM.nobs);
    OptionsMoM.variance_correction_factor = (1+1/OptionsMoM.simulation_multiple);
    % draw shocks for SMM
    smmstream = RandStream('mt19937ar','Seed',OptionsMoM.seed);
    temp_shocks = randn(smmstream,OptionsMoM.long+OptionsMoM.drop,Model.exo_nbr);
    if OptionsMoM.bounded_shock_support == 1
        temp_shocks(temp_shocks>2) = 2;
        temp_shocks(temp_shocks<-2) = -2;
    end
    OptionsMoM.shock_series = temp_shocks;
end

% -------------------------------------------------------------------------
% Step 5: checks for steady state at initial parameters
% -------------------------------------------------------------------------

% Check for logged steady state ...is this necessary @wmutschl
if OptionsMoM.logged_steady_state
    DynareResults.dr.ys=exp(DynareResults.dr.ys);
    DynareResults.steady_state=exp(DynareResults.steady_state);
    OptionsMoM.logged_steady_state=0;
end

% setting steadystate_check_flag option
if OptionsMoM.steadystate.nocheck
    steadystate_check_flag = 0;
else
    steadystate_check_flag = 1;
end

% Check steady state at initial model parameter values
[DynareResults.steady_state, params, info] = evaluate_steady_state(DynareResults.steady_state,Model,OptionsMoM,DynareResults,steadystate_check_flag);
if info(1)
    fprintf('\nmethod_of_moments: The steady state at the initial parameters cannot be computed.\n')
    print_info(info, 0, OptionsMoM);
end

try
    % check if steady state solves static model
    [DynareResults.steady_state, new_steady_params] = evaluate_steady_state(DynareResults.steady_state,Model,OptionsMoM,DynareResults,1);

    % check whether steady state file changes estimated parameters
    if isfield(EstimatedParameters,'param_vals') && ~isempty(EstimatedParameters.param_vals)
        Model0=Model; %store Model structure
        old_steady_params=Model.params; %save initial parameters for check if steady state changes param values

        Model0.params(EstimatedParameters.param_vals(:,1))=Model0.params(EstimatedParameters.param_vals(:,1))*1.01; %vary parameters
        [~, new_steady_params_2] = evaluate_steady_state(DynareResults.steady_state,Model0,OptionsMoM,DynareResults,1);

        changed_par_indices=find((old_steady_params(EstimatedParameters.param_vals(:,1))-new_steady_params(EstimatedParameters.param_vals(:,1))) ...
                                 | (Model0.params(EstimatedParameters.param_vals(:,1))-new_steady_params_2(EstimatedParameters.param_vals(:,1))));

        if ~isempty(changed_par_indices)
            fprintf('\nThe steady state file internally changed the values of the following estimated parameters:\n')
            disp(char(Model.param_names(EstimatedParameters.param_vals(changed_par_indices,1))))
            fprintf('This will override parameter values and may lead to wrong results.\n')
            fprintf('Check whether this is really intended.\n')
            warning('The steady state file internally changes the values of the estimated parameters.')
        end
    end

    % display warning if some parameters are still NaN
    test_for_deep_parameters_calibration(Model);

catch % if check fails, provide info on using calibration if present
    e = lasterror();
    if EstimatedParameters.full_calibration_detected %calibrated model present and no explicit starting values
        skipline(1);
        fprintf('There was an error in computing the moments for initial parameter values.\n')
        fprintf('If this is not a problem with the setting of options (check the error message below),\n')
        fprintf('you should try using the calibrated version of the model as starting values. To do\n')
        fprintf('this, add an empty estimated_params_init-block with use_calibration option immediately before the estimation\n')
        fprintf('command (and after the estimated_params-block so that it does not get overwritten):\n');
        skipline(2);
    end
    rethrow(e);
end

% If steady state of observed variables is non zero, set noconstant equal 0
if (~OptionsMoM.loglinear && all(abs(DynareResults.steady_state(DynareResults.dr.order_var(DynareResults.dr.obs_var)))<1e-9)) || (OptionsMoM.loglinear && all(abs(log(DynareResults.steady_state(DynareResults.dr.order_var(DynareResults.dr.obs_var))))<1e-9))
    OptionsMoM.noconstant = 1;
else
    OptionsMoM.noconstant = 0;
    % If the data are prefiltered then there must not be constants in the
    % measurement equation of the DSGE model
    if OptionsMoM.prefilter
        skipline()
        disp('You have specified the option "prefilter" to demean your data but the')
        disp('steady state of of the observed variables is non zero.')
        disp('Either change the measurement equations, by centering the observed')
        disp('variables in the model block, or drop the prefiltering.')
        error('method_of_moments: The option ''prefilter'' is inconsistent with the non-zero mean measurement equations in the model.')
    end
end


% -------------------------------------------------------------------------
% Step 6: checks for objective function at initial parameters
% -------------------------------------------------------------------------
try
    % Check for NaN or complex values of moment-distance-funtion evaluated
    % at initial parameters and identity weighting matrix    
    DynareResults.mom.Sw = eye(OptionsMoM.mom_nbr);
    tic_id = tic;    
    [fval, info, exit_flag, DynareResults, Model, OptionsMoM] = feval(objective_function, xparam0, Bounds, DynareResults, EstimatedParameters, MatchedMoments, Model, OptionsMoM);
    if OptionsMoM.mode_compute == 13
        fval = fval'*fval;
    end
    elapsed_time = toc(tic_id);    
    if isnan(fval)
        error('method_of_moments: The initial value of the target function is NaN')
    elseif imag(fval)
        error('method_of_moments: The initial value of the target function is complex')
    end
    if info(1) > 0
        disp('method_of_moments: Error in computing the target function for initial parameter values')
        print_info(info, OptionsMoM.noprint, OptionsMoM)
    end
    if OptionsMoM.prefilter==1
        if (~OptionsMoM.loglinear && any(abs(DynareResults.steady_state(DynareResults.dr.order_var(DynareResults.dr.obs_var)))>1e-9)) || (OptionsMoM.loglinear && any(abs(log(DynareResults.steady_state(DynareResults.dr.order_var(DynareResults.dr.obs_var))))>1e-9))
            disp(['You are trying to estimate a model with a non zero steady state for the observed endogenous'])
            disp(['variables using demeaned data!'])
            error('method_of_moments: You should change something in your mod file...')
        end
    end    
    fprintf('Initial value of the moment objective function with identity weighting matrix: %6.4f \n\n', fval);
    fprintf('Time required to compute target function once: %5.4f seconds \n', elapsed_time);
    
catch % if check fails, provide info on using calibration if present
    e = lasterror();
    if EstimatedParameters.full_calibration_detected %calibrated model present and no explicit starting values
        skipline(1);
        fprintf('There was an error in computing the moments for initial parameter values.\n')
        fprintf('If this is not a problem with the setting of options (check the error message below),\n')
        fprintf('you should try using the calibrated version of the model as starting values. To do\n')
        fprintf('this, add an empty estimated_params_init-block with use_calibration option immediately before the estimation\n')
        fprintf('command (and after the estimated_params-block so that it does not get overwritten):\n');
        skipline(2);
    end
    rethrow(e);
end

if OptionsMoM.mode_compute == 0 %We only report value of moments distance at initial value of the parameters
    fprintf('no minimization of moments distance due to ''mode_compute=0''\n')
    return
end

% -------------------------------------------------------------------------
% Step 7a: Method of moments estimation: print some info
% -------------------------------------------------------------------------
fprintf('\n---------------------------------------------------\n')
if strcmp(OptionsMoM.mom_method,'SMM')
    fprintf('Simulated method of moments with');
elseif strcmp(OptionsMoM.mom_method,'GMM')
    fprintf('General method of moments with');
end
if OptionsMoM.prefilter
    fprintf('\n  - centered moments (prefilter=1)');
else
    fprintf('\n  - uncentered moments (prefilter=0)');
end
if OptionsMoM.penalized_estimator
    fprintf('\n  - penalized estimation using priors');
end
if     OptionsMoM.mode_compute ==   1; fprintf('\n  - optimizer (mode_compute=1): fmincon');
elseif OptionsMoM.mode_compute ==   2; fprintf('\n  - optimizer (mode_compute=2): continuous simulated annealing');
elseif OptionsMoM.mode_compute ==   3; fprintf('\n  - optimizer (mode_compute=3): fminunc');
elseif OptionsMoM.mode_compute ==   4; fprintf('\n  - optimizer (mode_compute=4): csminwel');
elseif OptionsMoM.mode_compute ==   7; fprintf('\n  - optimizer (mode_compute=7): fminsearch');
elseif OptionsMoM.mode_compute ==   8; fprintf('\n  - optimizer (mode_compute=8): Dynare Nelder-Mead simplex');
elseif OptionsMoM.mode_compute ==   9; fprintf('\n  - optimizer (mode_compute=9): CMA-ES');
elseif OptionsMoM.mode_compute ==  10; fprintf('\n  - optimizer (mode_compute=10): simpsa');
elseif OptionsMoM.mode_compute ==  12; fprintf('\n  - optimizer (mode_compute=12): particleswarm');
elseif OptionsMoM.mode_compute == 101; fprintf('\n  - optimizer (mode_compute=101): SolveOpt');
elseif OptionsMoM.mode_compute == 102; fprintf('\n  - optimizer (mode_compute=102): simulannealbnd');
elseif OptionsMoM.mode_compute ==  13; fprintf('\n  - optimizer (mode_compute=13): lsqnonlin');
end
if OptionsMoM.silent_optimizer
    fprintf(' (silent)');
end
fprintf('\n  - perturbation order:        %d', OptionsMoM.order)
if OptionsMoM.order > 1 && OptionsMoM.pruning
fprintf(' (with pruning)')
end
fprintf('\n  - number of matched moments: %d', OptionsMoM.mom_nbr);
fprintf('\n  - number of parameters:      %d', length(xparam0));
% Check if enough moments for estimation
if OptionsMoM.mom_nbr < length(xparam0)
    fprintf('\n');
    error('method_of_moments: We must have at least as many moments as parameters for a method of moments estimation.')
end
fprintf('\n\n')

% -------------------------------------------------------------------------
% Step 7b: Method of moments estimation: First-stage
% -------------------------------------------------------------------------
fprintf('First-stage estimation\n');
switch lower(OptionsMoM.weighting_matrix)
    case 'identity_matrix'
        fprintf('  - identity weighting matrix\n');
        DynareResults.mom.Sw = eye(length(DynareResults.mom.dataMoments));
    case 'diagonal'
        %@wmutschl: better description in fprintf
        fprintf('  - diagonal weighting matrix: diagonal of Newey-West estimator with lag order %d\n', OptionsMoM.bartlett_kernel_lag);
        fprintf('                      and data moments as estimate of unconditional moments\n');
        Wopt = method_of_moments_optimal_weighting_matrix(DynareResults.mom.m_data, DynareResults.mom.dataMoments, OptionsMoM.bartlett_kernel_lag);
        DynareResults.mom.Sw = chol(diag(diag(Wopt)));
    case 'optimal'
        %@wmutschl: better description in fprintf
        fprintf('  - weighting matrix: optimal. At first-stage we use diagonal of Newey-West estimator with lag order %d\n', OptionsMoM.bartlett_kernel_lag);
        fprintf('                      and the data moments as initial estimate of unconditional moments\n');
        Wopt = method_of_moments_optimal_weighting_matrix(DynareResults.mom.m_data, DynareResults.mom.dataMoments, OptionsMoM.bartlett_kernel_lag);
        DynareResults.mom.Sw = chol(diag(diag(Wopt)));
    otherwise %user specified matrix in file
        fprintf('  - weighting matrix: user-specified\n');
        try
            load(OptionsMoM.weighting_matrix,'weighting_matrix')            
        catch
            error(['method_of_moments: No matrix named ''weighting_matrix'' could be found in ',OptionsMoM.weighting_matrix,'.mat'])
        end
        [nrow, ncol] = size(weighting_matrix);
        if ~isequal(nrow,ncol) && ~isequal(nrow,length(DynareResults.mom.dataMoments)) %check if square and right size
            error(['method_of_moments: weighting_matrix must be square and have ',num2str(length(DynareResults.mom.dataMoments)),' rows and columns'])
        end
        try %check for positive definiteness
            DynareResults.Sw = chol(weighting_matrix);
            hsd = sqrt(diag(weighting_matrix));
            inv(weighting_matrix./(hsd*hsd'))./(hsd*hsd');
        catch
            error('method_of_moments: Specified weighting_matrix is not positive definite')
        end
end
Woptflag = 0;
xparam1 = xparam0;
for istep1 = 1:2
    [xparam1, fval, exitflag, hessian_mat, OptionsMoM] = dynare_minimize_objective(objective_function, xparam1, OptionsMoM.mode_compute, OptionsMoM, [Bounds.lb Bounds.ub], BayesInfo.name, BayesInfo, [],...
                                                                                   Bounds, DynareResults, EstimatedParameters, MatchedMoments, Model, OptionsMoM);
    if OptionsMoM.mode_compute == 13
        fval = fval'*fval;
    end
    fprintf('\nIteration %d value of minimized moment''s distance target function: %f.\n',istep1,fval)
    if OptionsMoM.verbose
        DynareResults.mom=display_estimation_results_table(xparam1,NaN(size(xparam1)),Model,OptionsMoM,EstimatedParameters,BayesInfo,DynareResults.mom,DynareResults.prior.pnames,sprintf('%s (FIRST-STAGE ITERATION %d) verbose',OptionsMoM.mom_method,istep1),sprintf('verbose_%s_1st_stage_iter_%d',lower(OptionsMoM.mom_method),istep1));
    end
end

% Update Model and DynareResults (in particular DynareResults.mom.modelMoments)
Model = set_all_parameters(xparam1,EstimatedParameters,Model);
[fval, ~, ~, DynareResults, ~, ~] = feval(objective_function, xparam1, Bounds, DynareResults, EstimatedParameters, MatchedMoments, Model, OptionsMoM);
if OptionsMoM.mode_compute == 13
    fval = fval'*fval;
end

% Compute Standard errors
SE_1 = method_of_moments_standard_errors(xparam1, objective_function, Bounds, DynareResults, EstimatedParameters, MatchedMoments, Model, OptionsMoM, Woptflag);

% Store first-stage results in output structure
DynareResults.mom = display_estimation_results_table(xparam1,SE_1,Model,OptionsMoM,EstimatedParameters,BayesInfo,DynareResults.mom,DynareResults.prior.pnames,sprintf('%s (FIRST-STAGE)',OptionsMoM.mom_method),sprintf('%s_1st_stage',lower(OptionsMoM.mom_method)));

% -------------------------------------------------------------------------
% Step 7c: Method of moments estimation: Second-stage
% -------------------------------------------------------------------------
fprintf('Second-stage estimation\n');
switch lower(OptionsMoM.weighting_matrix)
    case 'identity_matrix'
        fprintf('  - weighting matrix: identity\n');
        DynareResults.mom.Sw = eye(length(DynareResults.mom.dataMoments));
    case 'diagonal'
        fprintf('  - weighting matrix: diagonal of Newey-West estimator with lag order %d\n', OptionsMoM.bartlett_kernel_lag);
        fprintf('                      and based on first-stage estimate of unconditional model moments\n');
        Wopt = method_of_moments_optimal_weighting_matrix(DynareResults.mom.m_data, DynareResults.mom.modelMoments, OptionsMoM.bartlett_kernel_lag);
        DynareResults.mom.Sw = chol(diag(diag(Wopt)));
    case 'optimal'
        fprintf('  - weighting matrix: Newey-West estimator with lag order %d\n', OptionsMoM.bartlett_kernel_lag);
        fprintf('                      and based on first-stage estimate of unconditional model moments\n');
        Wopt = method_of_moments_optimal_weighting_matrix(DynareResults.mom.m_data, DynareResults.mom.modelMoments, OptionsMoM.bartlett_kernel_lag);
        DynareResults.mom.Sw = chol(Wopt);
        Woptflag = 1;
        fprintf('                      rank of optimal weighting matrix: %d\n',rank(Wopt));
    otherwise %keep user specified matrix in file
        fprintf('  - weighting matrix: user-specified\n');
end

xparam2 = xparam1;
for istep2 = 1:2
    [xparam2, fval, exitflag, hessian_mat, OptionsMoM] = dynare_minimize_objective(objective_function, xparam2, OptionsMoM.mode_compute, OptionsMoM, [Bounds.lb Bounds.ub], BayesInfo.name, BayesInfo, [],...
                                                                                   Bounds, DynareResults, EstimatedParameters, MatchedMoments, Model, OptionsMoM);
    if OptionsMoM.mode_compute == 13
        fval = fval'*fval;
    end
    fprintf('\n  - iteration %d value of minimized moment''s distance target function: %f.\n',istep2,fval)
    if OptionsMoM.verbose
        DynareResults.mom=display_estimation_results_table(xparam2,NaN(size(xparam2)),Model,OptionsMoM,EstimatedParameters,BayesInfo,DynareResults.mom,DynareResults.prior.pnames,sprintf('%s (SECOND-STAGE ITERATION %d) verbose',OptionsMoM.mom_method,istep2),sprintf('verbose_%s_2nd_stage_iter_%d',lower(OptionsMoM.mom_method),istep2));
    end
end

% Update Model and DynareResults (in particular DynareResults.mom.modelMoments)
Model = set_all_parameters(xparam2,EstimatedParameters,Model);
[fval, ~, ~, DynareResults, ~, ~] = feval(objective_function, xparam2, Bounds, DynareResults, EstimatedParameters, MatchedMoments, Model, OptionsMoM);
if OptionsMoM.mode_compute == 13
    fval = fval'*fval;
end

% Compute Standard errors
SE_2 = method_of_moments_standard_errors(xparam2, objective_function, Bounds, DynareResults, EstimatedParameters, MatchedMoments, Model, OptionsMoM, Woptflag);

% Store second-stage results in output structure
DynareResults.mom = display_estimation_results_table(xparam2,SE_2,Model,OptionsMoM,EstimatedParameters,BayesInfo,DynareResults.mom,DynareResults.prior.pnames,sprintf('%s (SECOND-STAGE)',OptionsMoM.mom_method),sprintf('%s_2nd_stage',lower(OptionsMoM.mom_method)));

% Compute J statistic
if strcmp(OptionsMoM.mom_method,'SMM')    
    Variance_correction_factor = OptionsMoM.variance_correction_factor;
elseif strcmp(OptionsMoM.mom_method,'GMM')
    Variance_correction_factor=1;
end
DynareResults.mom.J_test.j_stat          = DynareDataset.nobs*Variance_correction_factor*fval;
DynareResults.mom.J_test.degrees_freedom = length(DynareResults.mom.modelMoments)-length(xparam2);
DynareResults.mom.J_test.p_val           = 1-chi2cdf(DynareResults.mom.J_test.j_stat, DynareResults.mom.J_test.degrees_freedom);
fprintf('\n p-value of J-test: %f\n',DynareResults.mom.J_test.p_val)

fprintf('\n==== Method of Moments Estimation Completed ====\n\n')

% -------------------------------------------------------------------------
% Step 8: Clean up
% -------------------------------------------------------------------------
% restore warnings
warning('on','MATLAB:singularMatrix');
