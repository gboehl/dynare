function options_mom_ = default_option_mom_values(options_mom_, options_, dname, doBayesianEstimation)
% function options_mom_ = default_option_mom_values(options_mom_, options_, dname, doBayesianEstimation)

% Returns structure containing the options for method_of_moments command

% options_mom_ is local and contains default and user-specified values for
% all settings needed for the method of moments estimation. Some options,
% though, are set by the preprocessor into options_ and we copy these over.
% The idea is to be independent of options_ and have full control of the
% estimation instead of possibly having to deal with options chosen somewhere
% else in the mod file.

% =========================================================================
% INPUTS
%  o options_mom_:           [structure] information about all (user-specified and updated) settings used in estimation (options_mom_)
%  o options_:               [structure] information on global options
%  o dname:                  [string]    name of directory to store results
%  o doBayesianEstimation    [boolean]   indicator whether we do Bayesian estimation
% -------------------------------------------------------------------------
% OUTPUTS
%  o oo_:                    [structure] storage for results (oo_)
%  o options_mom_:           [structure] information about all (user-specified and updated) settings used in estimation (options_mom_)
% -------------------------------------------------------------------------
% This function is called by
%  o mom.run
% -------------------------------------------------------------------------
% This function calls
% o set_default_option
% o user_has_matlab_license
% o user_has_octave_forge_package
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
% =========================================================================


mom_method = options_mom_.mom.mom_method; % this is a required option

% -------------------------------------------------------------------------
% LIMITATIONS
% -------------------------------------------------------------------------

if options_.logged_steady_state || options_.loglinear
    error('method_of_moments: The loglinear option is not supported. Please append the required logged variables as auxiliary equations.')
else
    options_mom_.logged_steady_state = 0;
    options_mom_.loglinear = false;
end
options_mom_.hessian.use_penalized_objective = false; % penalized objective not yet
% options related to variable declarations
if isfield(options_,'trend_coeffs')
    error('method_of_moments: %s does not allow for trend in data',mom_method)
end
% options related to endogenous prior restrictions are not supported
if ~isempty(options_.endogenous_prior_restrictions.irf) && ~isempty(options_.endogenous_prior_restrictions.moment)
    fprintf('method_of_moments: Endogenous prior restrictions are not supported yet and will be skipped.\n')
end
options_mom_.endogenous_prior_restrictions.irf    = {};
options_mom_.endogenous_prior_restrictions.moment = {};

options_mom_.mom.analytic_jacobian_optimizers = [1, 3, 4, 13, 101]; % these are currently supported optimizers that are able to use the analytical_jacobian option

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
%options_mom_ = set_default_option(options_mom_,'verbosity',false);     % 
if doBayesianEstimation
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
if (isoctave && user_has_octave_forge_package('optim')) || (~isoctave && user_has_matlab_license('optimization_toolbox'))
    if strcmp(mom_method,'GMM') || strcmp(mom_method,'SMM')
        options_mom_ = set_default_option(options_mom_,'mode_compute',13); % specifies lsqnonlin as default optimizer for minimization
    end
else
    options_mom_ = set_default_option(options_mom_,'mode_compute',4); % specifies csminwel as fallback default option for minimization
end
options_mom_ = set_default_option(options_mom_,'additional_optimizer_steps',[]);   % vector of additional mode-finders run after mode_compute
options_mom_ = set_default_option(options_mom_,'optim_opt',[]);                    % a list of NAME and VALUE pairs to set options for the optimization routines. Available options depend on mode_compute
options_mom_ = set_default_option(options_mom_,'silent_optimizer',false);          % run minimization of moments distance silently without displaying results or saving files in between
options_mom_ = set_default_option(options_mom_,'huge_number',1e7);                 % value for replacing the infinite bounds on parameters by finite numbers. Used by some optimizers for numerical reasons
options_mom_.mom = set_default_option(options_mom_.mom,'analytic_jacobian',false); % use analytic Jacobian in optimization, only available for GMM and gradient-based optimizers
options_mom_.optimizer_vec = [options_mom_.mode_compute;num2cell(options_mom_.additional_optimizer_steps)];

% perturbation related
options_mom_ = set_default_option(options_mom_,'order',1);                              % order of Taylor approximation in perturbation
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
options_mom_ = set_default_option(options_mom_,'sylvester_fp',false);                   % determines whether to use fixed point algorihtm to solve Sylvester equation (gensylv_fp), faster for large scale models
options_mom_ = set_default_option(options_mom_,'sylvester_fixed_point_tol',1e-12);      % convergence criterion used in the fixed point Sylvester solver

% mode check plot
options_mom_.mode_check.nolik = false;                                                          % we don't do likelihood (also this initializes mode_check substructure)
options_mom_.mode_check = set_default_option(options_mom_.mode_check,'status',false);           % plot the target function for values around the computed minimum for each estimated parameter in turn. This is helpful to diagnose problems with the optimizer.
options_mom_.mode_check = set_default_option(options_mom_.mode_check,'neighbourhood_size',.5);  % width of the window around the computed minimum to be displayed on the diagnostic plots. This width is expressed in percentage deviation. The Inf value is allowed, and will trigger a plot over the entire domain
options_mom_.mode_check = set_default_option(options_mom_.mode_check,'symmetric_plots',true);   % ensure that the check plots are symmetric around the minimum. A value of 0 allows to have asymmetric plots, which can be useful if the minimum is close to a domain boundary, or in conjunction with neighbourhood_size = Inf when the domain is not the entire real line
options_mom_.mode_check = set_default_option(options_mom_.mode_check,'number_of_points',20);    % number of points around the minimum where the target function is evaluated (for each parameter)


% -------------------------------------------------------------------------
% OPTIONS THAT NEED TO BE CARRIED OVER (E.G. SET BY THE PREPROCESSOR)
% -------------------------------------------------------------------------

% related to VAROBS block
options_mom_.varobs = options_.varobs;              % observable variables in order they are declared in varobs
options_mom_.varobs_id = options_.varobs_id;        % index for observable variables in M_.endo_names
options_mom_.obs_nbr = length(options_mom_.varobs); % number of observed variables

% related to call of dynare
options_mom_.console_mode = options_.console_mode;
options_mom_.parallel = options_.parallel;
options_mom_.parallel_info = options_.parallel_info;

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

% other
options_mom_.MaxNumberOfBytes = options_.MaxNumberOfBytes;
%options_mom_.MaximumNumberOfMegaBytes = options_.MaximumNumberOfMegaBytes;


% -------------------------------------------------------------------------
% DEFAULT VALUES
% -------------------------------------------------------------------------

options_mom_.analytic_derivation = 0;
options_mom_.analytic_derivation_mode = 0; % needed by get_perturbation_params_derivs.m, ie use efficient sylvester equation method to compute analytical derivatives as in Ratto & Iskrev (2012)
options_mom_.initialize_estimated_parameters_with_the_prior_mode = 0; % needed by set_prior.m
options_mom_.figures = options_.figures; % needed by plot_priors.m
options_mom_.ramsey_policy = false; % needed by evaluate_steady_state
options_mom_.risky_steadystate = false;  % needed by resol
options_mom_.jacobian_flag = true; % needed by dynare_solve