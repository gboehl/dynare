function [oo_, options_mom_, M_] = method_of_moments(bayestopt_, options_, oo_, estim_params_, M_, matched_moments_, options_mom_)
%function [oo_, options_mom_, M_] = method_of_moments(bayestopt_, options_, oo_, estim_params_, M_, matched_moments_, options_mom_)

% -------------------------------------------------------------------------
% This function performs a method of moments estimation with the following steps:
%   Step 0: Check if required structures and options exist
%   Step 1: - Prepare options_mom_ structure
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
%  o bayestopt_:             [structure] information about priors
%  o options_:               [structure] information about global options
%  o oo_:                    [structure] storage for results
%  o estim_params_:          [structure] information about estimated parameters
%  o M_:                     [structure] information about model
%  o matched_moments_:       [cell] information about selected moments to match in estimation
%                                         vars: matched_moments_{:,1});
%                                         lead/lags: matched_moments_{:,2}; 
%                                         powers: matched_moments_{:,3};
%  o options_mom_:           [structure] information about settings specified by the user
% -------------------------------------------------------------------------
% OUTPUTS
%  o oo_:                    [structure] storage for results (oo_)
%  o options_mom_:           [structure] information about all used settings used in estimation (options_mom_)
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
%  o makedataset.m
%  o method_of_moments_data_moments.m
%  o plot_priors.m
%  o print_info.m
%  o prior_bounds.m
%  o set_default_option.m
%  o set_prior.m
%  o set_state_space.m
%  o set_all_parameters.m
%  o test_for_deep_parameters_calibration.m
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
% - [ ] why does lsqnonlin take less time in Andreasen toolbox?
% - [ ] test user-specified weightning matrix
% - [ ] which qz_criterium value?
% - [ ] document that in method_of_moments_data_moments.m NaN are replaced by mean of moment
% - [ ] add IRF matching
% - [ ] test estimated_params_bounds block
% - [ ] test what happens if all parameters will be estimated but some/all are not calibrated
% - [ ] speed up lyapunov equation by using doubling with old initial values
% - [ ] check smm at order > 3 without pruning
% - [ ] provide option to use analytical derivatives to compute std errors (similar to what we already do in identification)
% - [ ] add Bayesian GMM/SMM estimation
% - [ ] useautocorr

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
if isempty(matched_moments_) % structure storing the moments used for the method of moments estimation
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

fprintf('\n==== Method of Moments (%s) Estimation ====\n\n',options_mom_.mom.mom_method)

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
    options_mom_.mom = set_default_option(options_mom_.mom,'bartlett_kernel_lag',20);             % bandwith in optimal weighting matrix
    options_mom_.mom = set_default_option(options_mom_.mom,'penalized_estimator',false);          % @wmutschl: provide description
    options_mom_.mom = set_default_option(options_mom_.mom,'verbose',false);                        % display and store intermediate estimation results
    options_mom_.mom = set_default_option(options_mom_.mom,'weighting_matrix','identity_matrix');   % weighting matrix in moments distance objective function; possible values: OPTIMAL|IDENTITY_MATRIX|DIAGONAL|filename
    options_mom_.mom = set_default_option(options_mom_.mom,'weighting_matrix_scaling_factor',1000); % scaling of weighting matrix
    options_mom_.mom = set_default_option(options_mom_.mom,'se_tolx',1e-5);           % step size for standard error
    options_mom_ = set_default_option(options_mom_,'order',1);                            % order of Taylor approximation in perturbation
    options_mom_ = set_default_option(options_mom_,'pruning',true);                       % use pruned state space system at higher-order
    % Checks for perturbation order
    if options_mom_.order < 1
        error('method_of_moments:: The order of the Taylor approximation cannot be 0!')
    end
end
if strcmp(options_mom_.mom.mom_method,'SMM')
    options_mom_.mom = set_default_option(options_mom_.mom,'burnin',500);                           % number of periods dropped at beginning of simulation
    options_mom_.mom = set_default_option(options_mom_.mom,'bounded_shock_support',false);        % trim shocks in simulation to +- 2 stdev
    options_mom_.mom = set_default_option(options_mom_.mom,'seed',24051986);                      % seed used in simulations
    options_mom_.mom = set_default_option(options_mom_.mom,'simulation_multiple',5);              % multiple of the data length used for simulation
    if options_mom_.mom.simulation_multiple < 1
        fprintf('The simulation horizon is shorter than the data. Dynare resets the simulation_multiple to 2.\n')
        options_mom_.mom.simulation_multiple = 2;
    end
end
if strcmp(options_mom_.mom.mom_method,'GMM')
    % Check for pruning with GMM at higher order
    if options_mom_.order > 1 && ~options_mom_.pruning
        fprintf('GMM at higher order only works with pruning, so we set pruning option to 1.\n');
        options_mom_.pruning = true;
    end
end
    
% General options that can be set by the user in the mod file, otherwise default values are provided
options_mom_ = set_default_option(options_mom_,'dirname',M_.fname); % directory in which to store estimation output
options_mom_ = set_default_option(options_mom_,'graph_format','eps');  % specify the file format(s) for graphs saved to disk
options_mom_ = set_default_option(options_mom_,'nodisplay',false);     % do not display the graphs, but still save them to disk
options_mom_ = set_default_option(options_mom_,'nograph',false);       % do not create graphs (which implies that they are not saved to the disk nor displayed)
options_mom_ = set_default_option(options_mom_,'noprint',false);       % do not print output to console
options_mom_ = set_default_option(options_mom_,'plot_priors',true);    % control plotting of priors
options_mom_ = set_default_option(options_mom_,'prior_trunc',1e-10);   % probability of extreme values of the prior density that is ignored when computing bounds for the parameters
options_mom_ = set_default_option(options_mom_,'TeX',false);           % print TeX tables and graphics

% Data and model options that can be set by the user in the mod file, otherwise default values are provided
options_mom_ = set_default_option(options_mom_,'first_obs',1);     % number of first observation
options_mom_ = set_default_option(options_mom_,'logdata',false);   % if data is already in logs
options_mom_ = set_default_option(options_mom_,'nobs',NaN);        % number of observations
options_mom_ = set_default_option(options_mom_,'prefilter',false); % demean each data series by its empirical mean and use centered moments
options_mom_ = set_default_option(options_mom_,'xls_sheet',1);     % name of sheet with data in Excel
options_mom_ = set_default_option(options_mom_,'xls_range','');    % range of data in Excel sheet
% Recursive estimation and forecast are not supported
if numel(options_mom_.nobs)>1
    error('method_of_moments: Recursive estimation and forecast for samples is not supported. Please set an integer as ''nobs''.');
end
if numel(options_mom_.first_obs)>1
    error('method_of_moments: Recursive estimation and forecast for samples is not supported. Please set an integer as ''first_obs''.');
end

% Optimization options that can be set by the user in the mod file, otherwise default values are provided
options_mom_ = set_default_option(options_mom_,'huge_number',1e7);               % value for replacing the infinite bounds on parameters by finite numbers. Used by some optimizers for numerical reasons
options_mom_ = set_default_option(options_mom_,'mode_compute',13);               % specifies the optimizer for minimization of moments distance
options_mom_ = set_default_option(options_mom_,'additional_optimizer_steps',[]); % vector of additional mode-finders run after mode_compute
options_mom_ = set_default_option(options_mom_,'optim_opt',[]);                  % a list of NAME and VALUE pairs to set options for the optimization routines. Available options depend on mode_compute
options_mom_ = set_default_option(options_mom_,'silent_optimizer',false);        % run minimization of moments distance silently without displaying results or saving files in between
% Mode_check plot options that can be set by the user in the mod file, otherwise default values are provided
options_mom_.mode_check.nolik = false;                                                          % we don't do likelihood (also this initializes mode_check substructure)
options_mom_.mode_check = set_default_option(options_mom_.mode_check,'status',true);            % plot the target function for values around the computed mode for each estimated parameter in turn. This is helpful to diagnose problems with the optimizer.
options_mom_.mode_check = set_default_option(options_mom_.mode_check,'neighbourhood_size',.5);  % width of the window around the mode to be displayed on the diagnostic plots. This width is expressed in percentage deviation. The Inf value is allowed, and will trigger a plot over the entire domain
options_mom_.mode_check = set_default_option(options_mom_.mode_check,'symmetric_plots',true);   % ensure that the check plots are symmetric around the mode. A value of 0 allows to have asymmetric plots, which can be useful if the posterior mode is close to a domain boundary, or in conjunction with mode_check_neighbourhood_size = Inf when the domain is not the entire real line
options_mom_.mode_check = set_default_option(options_mom_.mode_check,'number_of_points',20);    % number of points around the mode where the target function is evaluated (for each parameter)

% Numerical algorithms options that can be set by the user in the mod file, otherwise default values are provided
options_mom_ = set_default_option(options_mom_,'aim_solver',false);                     % Use AIM algorithm to compute perturbation approximation
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
options_mom_ = set_default_option(options_mom_,'qz_criterium',1-1e-6);                  % value used to split stable from unstable eigenvalues in reordering the Generalized Schur decomposition used for solving first order problems [IS THIS CORRET @wmutschl]
options_mom_ = set_default_option(options_mom_,'qz_zero_threshold',1e-6);               % value used to test if a generalized eigenvalue is 0/0 in the generalized Schur decomposition
if options_mom_.order > 2
    fprintf('Dynare will use ''k_order_solver'' as the order>2\n');
    options_mom_.k_order_solver = true;
end

% -------------------------------------------------------------------------
% Step 1b: Options that are set by the preprocessor and (probably) need to be carried over
% -------------------------------------------------------------------------

% options related to VAROBS
if ~isfield(options_,'varobs')
    error('method_of_moments: VAROBS statement is missing!')
else
    options_mom_.varobs  = options_.varobs;     % observable variables in declaration order
    options_mom_.obs_nbr = length(options_mom_.varobs);% number of observed variables
    % Check that each declared observed variable is also an endogenous variable
    for i = 1:options_mom_.obs_nbr
        if ~any(strcmp(options_mom_.varobs{i},M_.endo_names))
            error(['method_of_moments: Unknown variable (' options_mom_.varobs{i} ')!'])
        end
    end

    % Check that a variable is not declared as observed more than once.
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

% options related to model; block
options_mom_.linear   = options_.linear;
options_mom_.use_dll  = options_.use_dll;
options_mom_.block    = options_.block;
options_mom_.bytecode = options_.bytecode;

% options related to steady; command
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

% options related to endogenous prior restrictions
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

options_mom_.vector_output= false;           % specifies whether the objective function returns a vector

% -------------------------------------------------------------------------
% Step 1d: Other options that need to be initialized
% -------------------------------------------------------------------------
options_mom_.initialize_estimated_parameters_with_the_prior_mode = 0; % needed by set_prior.m
options_mom_.figures.textwidth = 0.8; %needed by plot_priors.m
options_mom_.ramsey_policy = 0; % needed by evaluate_steady_state
options_mom_.debug = false; %needed by resol.m
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
% Step 2: Checks and transformations for matched moments structure (preliminary)
% -------------------------------------------------------------------------
% Note that we do not have a preprocessor interface yet for this, so this
% will need much improvement later on. @wmutschl

% Initialize indices
options_mom_.mom.index.E_y       = false(options_mom_.obs_nbr,1);                      %unconditional first order product moments
options_mom_.mom.index.E_yy      = false(options_mom_.obs_nbr,options_mom_.obs_nbr);   %unconditional second order product moments
options_mom_.mom.index.E_yyt     = false(options_mom_.obs_nbr,options_mom_.obs_nbr,0); %unconditional temporal second order product moments
options_mom_.mom.index.E_y_pos   = zeros(options_mom_.obs_nbr,1);                      %position in matched moments block
options_mom_.mom.index.E_yy_pos  = zeros(options_mom_.obs_nbr,options_mom_.obs_nbr);   %position in matched moments block
options_mom_.mom.index.E_yyt_pos = zeros(options_mom_.obs_nbr,options_mom_.obs_nbr,0); %position in matched moments block

for jm=1:size(matched_moments_,1)
    % higher-order product moments not supported yet for GMM
    if strcmp(options_mom_.mom.mom_method, 'GMM') && sum(matched_moments_{jm,3}) > 2
        error('method_of_moments: GMM does not yet support product moments higher than 2. Change row %d in ''matched_moments'' block.',jm);
    end    
    % Check if declared variables are also observed (needed as otherwise the dataset variables won't coincide)
    if any(~ismember(oo_.dr.inv_order_var(matched_moments_{jm,1})', oo_.dr.obs_var))
        error('method_of_moments: Variables in row %d in ''matched_moments'' block need to be declared as VAROBS.', jm)
    end
    
    if strcmp(options_mom_.mom.mom_method, 'GMM')
    % Check (for now) that only lags are declared
        if any(matched_moments_{jm,2}>0)
            error('method_of_moments: Leads in row %d in the ''matched_moments'' block are not supported for GMM, shift the moments and declare only lags.', jm)
        end
        % Check (for now) that first declared variable has zero lag
        if matched_moments_{jm,2}(1)~=0
            error('method_of_moments: The first variable declared in row %d in the ''matched_moments'' block is not allowed to have a lead or lag for GMM;\n                   reorder the variables in the row such that the first variable has zero lag!',jm)
        end
    end
    vars = oo_.dr.inv_order_var(matched_moments_{jm,1})';
    if sum(matched_moments_{jm,3}) == 1
        % First-order product moment
        vpos = (oo_.dr.obs_var == vars);
        options_mom_.mom.index.E_y(vpos,1) = true;
        options_mom_.mom.index.E_y_pos(vpos,1) = jm;
        matched_moments_{jm,4}=['E(',M_.endo_names{matched_moments_{jm,1}},')'];
        matched_moments_{jm,5}=['$E(',M_.endo_names_tex{matched_moments_{jm,1}},')$'];
    elseif sum(matched_moments_{jm,3}) == 2
        % Second-order product moment
        idx1 = (oo_.dr.obs_var == vars(1));
        idx2 = (oo_.dr.obs_var == vars(2));
        lag1 = matched_moments_{jm,2}(1);
        lag2 = matched_moments_{jm,2}(2);
        if lag1==0 && lag2==0 % contemporaneous covariance matrix
            options_mom_.mom.index.E_yy(idx1,idx2) = true;
            options_mom_.mom.index.E_yy(idx2,idx1) = true;
            options_mom_.mom.index.E_yy_pos(idx1,idx2) = jm;
            options_mom_.mom.index.E_yy_pos(idx2,idx1) = jm;
            matched_moments_{jm,4}=['E(',M_.endo_names{matched_moments_{jm,1}(1)},',',M_.endo_names{matched_moments_{jm,1}(2)},')'];
            matched_moments_{jm,5}=['$E({',M_.endo_names_tex{matched_moments_{jm,1}(1)},'}_t,{',M_.endo_names_tex{matched_moments_{jm,1}(1)},'}_t)$'];
        elseif lag1==0 && lag2 < 0
            options_mom_.mom.index.E_yyt(idx1,idx2,-lag2) = true;
            options_mom_.mom.index.E_yyt_pos(idx1,idx2,-lag2) = jm;
            matched_moments_{jm,4}=['E(',M_.endo_names{matched_moments_{jm,1}(1)},',',M_.endo_names{matched_moments_{jm,1}(2)},'(',num2str(lag2),'))'];
            matched_moments_{jm,5}=['$E({',M_.endo_names_tex{matched_moments_{jm,1}(1)},'}_t\times{',M_.endo_names_tex{matched_moments_{jm,1}(1)},'_{t',num2str(lag2) ,'})$'];
        end
    end
end


% @wmutschl: add check for duplicate moments by using the cellfun and unique functions
%Remove duplicate elements
UniqueMomIdx = [nonzeros(options_mom_.mom.index.E_y_pos); nonzeros(tril(options_mom_.mom.index.E_yy_pos)); nonzeros(options_mom_.mom.index.E_yyt_pos)];
DuplicateMoms = setdiff(1:size(matched_moments_,1),UniqueMomIdx);
if ~isempty(DuplicateMoms)
    fprintf('Found and removed duplicate declared moments in ''matched_moments'' block in rows: %s.\n',num2str(DuplicateMoms))
end
%reorder matched_moments_ to be compatible with options_mom_.mom.index
matched_moments_ = matched_moments_(UniqueMomIdx,:);
if strcmp(options_mom_.mom.mom_method,'SMM')
    options_mom_.mom=rmfield(options_mom_.mom,'index');
end

% Check if both prefilter and first moments were specified
options_mom_.mom.first_moment_indicator = find(cellfun(@(x) sum(abs(x))==1,matched_moments_(:,3)))';
if options_mom_.prefilter && ~isempty(options_mom_.mom.first_moment_indicator)
    fprintf('Centered moments requested (prefilter option is set); therefore, ignore declared first moments in ''matched_moments'' block in rows: %u.\n',options_mom_.mom.first_moment_indicator');
    matched_moments_(options_mom_.mom.first_moment_indicator,:)=[]; %remove first moments entries
    options_mom_.mom.first_moment_indicator = [];
end
options_mom_.mom.mom_nbr = size(matched_moments_,1);

% Get maximum lag number for autocovariances/autocorrelations
options_mom_.ar = max(cellfun(@max,matched_moments_(:,2))) - min(cellfun(@min,matched_moments_(:,2)));

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
if any(bayestopt_laplace.pshape > 0) % prior specified, not ML
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
        warning('off','MATLAB:singularMatrix');
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

% Check whether on initialization
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

% missing observations are not supported yet
if any(any(isnan(dataset_.data)))
    error('method_of_moments: missing observations are not supported')
end

% Check length of data for estimation of second moments
if options_mom_.ar > options_mom_.nobs+1
    error('method_of_moments: Data set is too short to compute second moments');
end

% Get data moments for the method of moments
[oo_.mom.data_moments, oo_.mom.m_data] = method_of_moments_data_moments(dataset_.data, oo_, matched_moments_, options_mom_);

% Get shock series for SMM and set variance correction factor
if strcmp(options_mom_.mom.mom_method,'SMM')
    options_mom_.mom.long = round(options_mom_.mom.simulation_multiple*options_mom_.nobs);
    options_mom_.mom.variance_correction_factor = (1+1/options_mom_.mom.simulation_multiple);
    % draw shocks for SMM
    smmstream = RandStream('mt19937ar','Seed',options_mom_.mom.seed);
    temp_shocks = randn(smmstream,options_mom_.mom.long+options_mom_.mom.burnin,M_.exo_nbr);
    temp_shocks_ME = randn(smmstream,options_mom_.mom.long,length(M_.H));
    if options_mom_.mom.bounded_shock_support == 1
        temp_shocks(temp_shocks>2) = 2;
        temp_shocks(temp_shocks<-2) = -2;
        temp_shocks_ME(temp_shocks_ME<-2) = -2;
        temp_shocks_ME(temp_shocks_ME<-2) = -2;
    end
    options_mom_.mom.shock_series = temp_shocks;
    options_mom_.mom.ME_shock_series = temp_shocks_ME;
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
    end
end

% display warning if some parameters are still NaN
test_for_deep_parameters_calibration(M_);

% If steady state of observed variables is non zero, set noconstant equal 0
if all(abs(oo_.steady_state(oo_.dr.order_var(oo_.dr.obs_var)))<1e-9)
    options_mom_.noconstant = 1;
else
    options_mom_.noconstant = 0;
end

% -------------------------------------------------------------------------
% Step 6: checks for objective function at initial parameters
% -------------------------------------------------------------------------
objective_function = str2func('method_of_moments_objective_function');
try
    % Check for NaN or complex values of moment-distance-funtion evaluated
    % at initial parameters and identity weighting matrix    
    oo_.mom.Sw = eye(options_mom_.mom.mom_nbr);
    tic_id = tic;    
    [fval, info, ~, ~, ~, oo_, M_] = feval(objective_function, xparam0, Bounds, oo_, estim_params_, matched_moments_, M_, options_mom_);
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

if options_mom_.mode_compute == 0 %We only report value of moments distance at initial value of the parameters
    fprintf('No minimization of moments distance due to ''mode_compute=0''\n')
    return
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
    fprintf('\n  - penalized estimation using priors');
end
if     options_mom_.mode_compute ==   1; fprintf('\n  - optimizer (mode_compute=1): fmincon');
elseif options_mom_.mode_compute ==   2; fprintf('\n  - optimizer (mode_compute=2): continuous simulated annealing');
elseif options_mom_.mode_compute ==   3; fprintf('\n  - optimizer (mode_compute=3): fminunc');
elseif options_mom_.mode_compute ==   4; fprintf('\n  - optimizer (mode_compute=4): csminwel');
elseif options_mom_.mode_compute ==   5; fprintf('\n  - optimizer (mode_compute=5): newrat');
elseif options_mom_.mode_compute ==   6; fprintf('\n  - optimizer (mode_compute=6): gmhmaxlik');
elseif options_mom_.mode_compute ==   7; fprintf('\n  - optimizer (mode_compute=7): fminsearch');
elseif options_mom_.mode_compute ==   8; fprintf('\n  - optimizer (mode_compute=8): Dynare Nelder-Mead simplex');
elseif options_mom_.mode_compute ==   9; fprintf('\n  - optimizer (mode_compute=9): CMA-ES');
elseif options_mom_.mode_compute ==  10; fprintf('\n  - optimizer (mode_compute=10): simpsa');
elseif options_mom_.mode_compute ==  11; fprintf('\n  - optimizer (mode_compute=11): online_auxiliary_filter');
elseif options_mom_.mode_compute ==  12; fprintf('\n  - optimizer (mode_compute=12): particleswarm');
elseif options_mom_.mode_compute == 101; fprintf('\n  - optimizer (mode_compute=101): SolveOpt');
elseif options_mom_.mode_compute == 102; fprintf('\n  - optimizer (mode_compute=102): simulannealbnd');
elseif options_mom_.mode_compute ==  13; fprintf('\n  - optimizer (mode_compute=13): lsqnonlin');
elseif ischar(minimizer_algorithm); fprintf(['\n  - user-defined optimizer: ' minimizer_algorithm]);
else
    error('method_of_moments: Unknown optimizer, please contact the developers ')
end
if options_mom_.silent_optimizer
    fprintf(' (silent)');
end
fprintf('\n  - perturbation order:        %d', options_mom_.order)
if options_mom_.order > 1 && options_mom_.pruning
    fprintf(' (with pruning)')
end
fprintf('\n  - number of matched moments: %d', options_mom_.mom.mom_nbr);
fprintf('\n  - number of parameters:      %d\n\n', length(xparam0));

% -------------------------------------------------------------------------
% Step 7b: Method of moments estimation: First-stage
% -------------------------------------------------------------------------
if size(options_mom_.mom.weighting_matrix,1)>1 && ~(any(strcmpi('diagonal',options_mom_.mom.weighting_matrix)) || any(strcmpi('optimal',options_mom_.mom.weighting_matrix)))
    fprintf('\nYou did not specify the use of an optimal or diagonal covariance matrix. There is no point in running an iterated method of moments.\n')
end
for stage_iter=1:size(options_mom_.mom.weighting_matrix,1)    
    Woptflag = 0;
    fprintf('Estimation stage %u\n',stage_iter);
    switch lower(options_mom_.mom.weighting_matrix{stage_iter})
        case 'identity_matrix'
            fprintf('  - identity weighting matrix\n');
            oo_.mom.Sw = eye(options_mom_.mom_nbr);
        case 'diagonal'
            %@wmutschl: better description in fprintf
            fprintf('  - diagonal weighting matrix: diagonal of Newey-West estimator with lag order %d\n', options_mom_.mom.bartlett_kernel_lag);
            fprintf('                      and data moments as estimate of unconditional moments\n');
            W_opt = method_of_moments_optimal_weighting_matrix(oo_.mom.m_data, oo_.mom.data_moments, options_mom_.mom.bartlett_kernel_lag);
            oo_.mom.Sw = chol(diag(diag(W_opt)));
        case 'optimal'
            %@wmutschl: better description in fprintf
            fprintf('  - weighting matrix: optimal. At first-stage we use diagonal of Newey-West estimator with lag order %d\n', options_mom_.mom.bartlett_kernel_lag);
            fprintf('                      and the data moments as initial estimate of unconditional moments\n');
            W_opt = method_of_moments_optimal_weighting_matrix(oo_.mom.m_data, oo_.mom.data_moments, options_mom_.mom.bartlett_kernel_lag);
            oo_.mom.Sw = chol(diag(diag(W_opt)));
            Woptflag = 1;
        otherwise %user specified matrix in file
            fprintf('  - weighting matrix: user-specified\n');
            try
                load(options_mom_.mom.weighting_matrix{stage_iter},'weighting_matrix')
            catch
                error(['method_of_moments: No matrix named ''weighting_matrix'' could be found in ',options_mom_.mom.weighting_matrix{stage_iter},'.mat'])
            end
            [nrow, ncol] = size(weighting_matrix);
            if ~isequal(nrow,ncol) && ~isequal(nrow,length(oo_.mom.data_moments)) %check if square and right size
                error(['method_of_moments: weighting_matrix must be square and have ',num2str(length(oo_.mom.data_moments)),' rows and columns'])
            end
            try %check for positive definiteness
                oo_.Sw = chol(weighting_matrix);
            catch
                error('method_of_moments: Specified weighting_matrix is not positive definite')
            end
    end
    
    optimizer_vec=[options_mom_.mode_compute,options_mom_.additional_optimizer_steps];
    
    for istep= 1:length(optimizer_vec)
        if optimizer_vec(istep)==13
            options_mom_.vector_output = true;
        else
            options_mom_.vector_output = false;
        end
        [xparam1, fval, exitflag] = dynare_minimize_objective(objective_function, xparam0, optimizer_vec(istep), options_mom_, [Bounds.lb Bounds.ub], bayestopt_laplace.name, bayestopt_laplace, [],...
            Bounds, oo_, estim_params_, matched_moments_, M_, options_mom_);
        if options_mom_.vector_output
            fval = fval'*fval;
        end
        fprintf('\nStage %d Iteration %d: value of minimized moment distance objective function: %12.10f.\n',stage_iter,istep,fval)
        if options_mom_.mom.verbose
            oo_.mom=display_estimation_results_table(xparam1,NaN(size(xparam1)),M_,options_mom_,estim_params_,bayestopt_laplace,oo_.mom,prior_dist_names,sprintf('%s (FIRST-STAGE ITERATION %d) verbose',options_mom_.mom.mom_method,istep),sprintf('verbose_%s_1st_stage_iter_%d',lower(options_mom_.mom.mom_method),istep));
        end
        xparam0=xparam1;
    end
    options_mom_.vector_output = false;
    
    % Update M_ and DynareResults (in particular oo_.mom.model_moments)
    M_ = set_all_parameters(xparam1,estim_params_,M_);
    [fval, ~, ~,~,~, oo_] = feval(objective_function, xparam1, Bounds, oo_, estim_params_, matched_moments_, M_, options_mom_);
    
    % Compute Standard errors
    SE = method_of_moments_standard_errors(xparam1, objective_function, Bounds, oo_, estim_params_, matched_moments_, M_, options_mom_, Woptflag);
    
    % Store first-stage results in output structure
    oo_.mom = display_estimation_results_table(xparam1,SE,M_,options_mom_,estim_params_,bayestopt_laplace,oo_.mom,prior_dist_names,sprintf('%s (STAGE %u)',options_mom_.mom.mom_method,stage_iter),sprintf('%s_stage_%u',lower(options_mom_.mom.mom_method),stage_iter));


end

%get optimal weighting matrix for J test, if necessary
if ~Woptflag
    W_opt = method_of_moments_optimal_weighting_matrix(oo_.mom.m_data, oo_.mom.model_moments, options_mom_.mom.bartlett_kernel_lag);
    oo_j=oo_;
    oo_j.mom.Sw = chol(W_opt);
    [fval] = feval(objective_function, xparam1, Bounds, oo_j, estim_params_, matched_moments_, M_, options_mom_);
end

% -------------------------------------------------------------------------
% Step 8: J test
% -------------------------------------------------------------------------
if options_mom_.mom.mom_nbr > length(xparam1)
    %get optimal weighting matrix for J test, if necessary
    if ~Woptflag
        W_opt = method_of_moments_optimal_weighting_matrix(oo_.mom.m_data, oo_.mom.model_moments, options_mom_.mom.bartlett_kernel_lag);
        oo_j=oo_;
        oo_j.mom.Sw = chol(W_opt);
        [fval] = feval(objective_function, xparam1, Bounds, oo_j, estim_params_, matched_moments_, M_, options_mom_);
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

title = ['Data moments and model moments (',options_mom_.mom.mom_method,')'];
headers = {'Moment','Data','Model','% dev. target'};
labels= matched_moments_(:,4);
data_mat=[oo_.mom.data_moments oo_.mom.model_moments 100*abs((oo_.mom.model_moments-oo_.mom.data_moments)./oo_.mom.data_moments)];
dyntable(options_mom_, title, headers, labels, data_mat, cellofchararraymaxlength(labels)+2, 10, 7);
if options_mom_.TeX
    lh = cellofchararraymaxlength(labels)+2;
    labels_TeX = matched_moments_(:,5);
    dyn_latex_table(M_, options_mom_, title, 'sim_corr_matrix', headers, labels_TeX, data_mat, lh, 10, 7);
end

if options_mom_.mode_check.status
    method_of_moments_mode_check(objective_function,xparam1,SE,options_mom_,M_,estim_params_,Bounds,bayestopt_laplace,...
        Bounds, oo_, estim_params_, matched_moments_, M_, options_mom_)
end

fprintf('\n==== Method of Moments Estimation (%s) Completed ====\n\n',options_mom_.mom.mom_method)

% -------------------------------------------------------------------------
% Step 8: Clean up
% -------------------------------------------------------------------------
% restore warnings
warning('on','MATLAB:singularMatrix');
