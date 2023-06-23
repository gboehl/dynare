function options_ = default_option_values(M_)
%function default_option_values()
% Returns structure containing the options for Dynare commands and their
% default values
%
% INPUTS
%    M_         [structure]     Model definition
%
% OUTPUTS
%    options    [structure]     Command options
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2018-2023 Dynare Team
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

options_.datafile = '';
options_.dirname = M_.fname;
options_.dataset = [];
options_.verbosity = 1;
options_.rplottype = 0;
options_.smpl = 0;
options_.dynatol.f = 1e-5;
options_.dynatol.x = 1e-5;
options_.timing = 0;
options_.gstep = ones(2,1);
options_.gstep(1) = 1e-2;
options_.gstep(2) = 1.0;
options_.jacobian_tolerance = []; %tolerance for rank of Jacobian in model_diagnostics
options_.debug = false;
options_.initval_file = false;
options_.schur_vec_tol = 1e-11; % used to find nonstationary variables in Schur decomposition of the
                                % transition matrix
options_.qz_criterium = [];
options_.qz_zero_threshold = 1e-6;
options_.lyapunov_complex_threshold = 1e-15;
options_.solve_algo = 4;
options_.solve_tolf = eps^(1/3);
options_.solve_tolx = eps^(2/3);
options_.trust_region_initial_step_bound_factor = 1;
options_.dr_display_tol=1e-6;
options_.minimal_workspace = false;
options_.dp.maxit = 3000;
options_.steady.maxit = 50;
options_.simul.maxit = 50;
options_.simul.robust_lin_solve = false;

options_.mode_check.status = false;
options_.mode_check.neighbourhood_size = .5;
options_.mode_check.symmetric_plots = true;
options_.mode_check.number_of_points = 20;
options_.mode_check.nolik = false;

options_.huge_number = 1e7;

% Default number of threads for parallelized mex files.
options_.threads.kronecker.sparse_hessian_times_B_kronecker_C = num_procs;
options_.threads.local_state_space_iteration_2 = num_procs;
options_.threads.local_state_space_iteration_3 = num_procs;
options_.threads.local_state_space_iteration_k = 1;
options_.threads.perfect_foresight_problem = num_procs;
options_.threads.k_order_perturbation = max(1, num_procs/2);

% steady state
options_.jacobian_flag = true;

% steady state file
if exist(['+' M_.fname '/steadystate.m'],'file')
    options_.steadystate_flag = 2;
elseif exist([M_.fname '_steadystate.m'],'file')
    options_.steadystate_flag = 1;
else
    options_.steadystate_flag = 0;
end
options_.steadystate_partial = [];
options_.steadystate.nocheck = false;

% subset of the estimated deep parameters
options_.ParamSubSet = 'None';

% bvar-dsge
options_.dsge_var = 0;
options_.dsge_varlag = 4;

% BVAR a la Sims
options_.bvar_replic = 2000;
options_.bvar_prior_tau = 3;
options_.bvar_prior_decay = 0.5;
options_.bvar_prior_lambda = 5;
options_.bvar_prior_mu = 2;
options_.bvar_prior_omega = 1;
options_.bvar_prior_flat = false;
options_.bvar_prior_train = 0;
options_.bvar.conf_sig = 0.6;

% Initialize the field that will contain the optimization algorithm's options declared in the
% estimation command (if any).
options_.optim_opt = [];

% Same for options to fsolve
options_.fsolve_options = [];

% Optimization algorithm [6] gmhmaxlik
gmhmaxlik.iterations = 3;
gmhmaxlik.number = 20000;
gmhmaxlik.nclimb = 200000;
gmhmaxlik.nscale = 200000;
gmhmaxlik.target = 1/3; % Target for the acceptance rate.
options_.gmhmaxlik = gmhmaxlik;

% Request user input.
options_.nointeractive = false;

% Graphics
options_.graphics.nrows = 3;
options_.graphics.ncols = 3;
options_.graphics.line_types = {'b-'};
options_.graphics.line_width = 1;
options_.graph_format = 'eps';
options_.nodisplay = false;
options_.nograph = false;
options_.no_graph.posterior = false;
options_.no_graph.shock_decomposition = false;
options_.no_graph.plot_shock_decomposition = false;
options_.XTick = [];
options_.XTickLabel = [];
options_.console_mode = false;
if isoctave
    if sum(get(0,'screensize'))==4
        options_.console_mode = true;
        options_.nodisplay = true;
    end
else
    if isunix && (~usejava('jvm') || ~feature('ShowFigureWindows'))
        options_.console_mode = true;
        options_.nodisplay = true;
    end
end

% IRFs & other stoch_simul output
options_.irf = 40;
options_.impulse_responses.plot_threshold=1e-10;
options_.zero_moments_tolerance=1e-10;
options_.relative_irf = false;
options_.ar = 5;
options_.hp_filter = 0;
options_.one_sided_hp_filter = 0;
options_.filtered_theoretical_moments_grid = 512;
options_.nodecomposition = false;
options_.nomoments = false;
options_.nocorr = false;
options_.periods = 0;
options_.noprint = false;
options_.SpectralDensity.trigger = false;
options_.SpectralDensity.plot  = 1;
options_.SpectralDensity.cutoff  = 150;
options_.SpectralDensity.sdl = 0.01;
options_.nofunctions = false;

options_.bandpass.indicator = false;
options_.bandpass.passband = [6; 32];
options_.bandpass.K=12;

options_.irf_opt.diagonal_only = false;
options_.irf_opt.stderr_multiples = false;
options_.irf_opt.irf_shock_graphtitles = {};
options_.irf_opt.irf_shocks = [];

% Extended path options
%
% Set debug flag
ep.debug = 0;
% Set memory flag
ep.memory = 0;
% Set verbose mode
ep.verbosity = 0;
% Set bytecode flag
ep.use_bytecode = 0;
% Initialization of the perfect foresight equilibrium paths
% * init=0, previous solution is used.
% * init=1, a path generated with the first order reduced form is used.
% * init=2, mix of cases 0 and 1.
ep.init = 0;
% Maximum number of iterations for the deterministic solver.
ep.maxit = 500;
% Number of periods for the perfect foresight model.
ep.periods = 200;
% Default step for increasing the number of periods if needed
ep.step = 50;
% Set check_stability flag
ep.check_stability = 0;
% Define last periods used to test if the solution is stable with respect to an increase in the number of periods.
ep.lp = 5;
% Define first periods used to test if the solution is stable with respect to an increase in the number of periods.
ep.fp = 2;
% Define the distribution for the structural innovations.
ep.innovation_distribution = 'gaussian';
% Set flag for the seed
ep.set_dynare_seed_to_default = 1;
% Set algorithm for the perfect foresight solver
ep.stack_solve_algo = 7;
ep.solve_algo = 9;
% Number of replications
ep.replic_nbr = 1;
% Parallel execution of replications
ep.parallel = false;
% Stochastic extended path related options.
ep.stochastic.IntegrationAlgorithm = 'Tensor-Gaussian-Quadrature'; % Other possible values are 'Stroud-Cubature-3' and 'Stroud-Cubature-5'
ep.stochastic.method = '';
ep.stochastic.algo = 0;
ep.stochastic.quadrature.ortpol = 'hermite';
ep.stochastic.order = 0;
ep.stochastic.quadrature.nodes = 5;
ep.stochastic.quadrature.pruned.status = 0;
ep.stochastic.quadrature.pruned.relative = 1e-5;
ep.stochastic.quadrature.pruned.level = 1e-5;
ep.stochastic.hybrid_order = 0;
% homotopic step in extended path simulations
ep.stochastic.homotopic_steps = true;
% Copy ep structure in options_ global structure
options_.ep = ep;


% Simulations of backward looking models options
%
bnlms.set_dynare_seed_to_default = true;
bnlms.innovation_distribution = 'gaussian';
options_.bnlms = bnlms;


% Particle filter
%
% Default is that we do not use the non linear kalman filter
particle.status = false;
% How do we initialize the states?
particle.initialization = 1;
particle.initial_state_prior_std = .1;
% Set the default number of particles.
particle.number_of_particles = 5000;
% Set the default approximation order (Smolyak)
particle.smolyak_accuracy = 3;
% By default we don't use pruning
particle.pruning = 0;
% Set the Gaussian approximation method for particles distributions
particle.approximation_method = 'unscented';
% Set unscented parameters alpha, beta and kappa for gaussian approximation
particle.unscented.alpha = 1;
particle.unscented.beta = 2;
particle.unscented.kappa = 1;
% Configuration of resampling in case of particles
particle.resampling.status.systematic = true;
particle.resampling.status.none = false;
particle.resampling.status.generic = false;
particle.resampling.threshold = .5;
particle.resampling.method.kitagawa = true;
particle.resampling.method.smooth = false;
particle.resampling.method.stratified = false;
% Set default algorithm
particle.filter_algorithm = 'sis';
% Approximation of the proposal distribution
particle.proposal_approximation.cubature = false;
particle.proposal_approximation.unscented = true;
particle.proposal_approximation.montecarlo = false;
% Approximation of the particle distribution
particle.distribution_approximation.cubature = false;
particle.distribution_approximation.unscented = true;
particle.distribution_approximation.montecarlo = false;
% Number of partitions for the smoothed resampling method
particle.resampling.number_of_partitions = 200;
% Configuration of the mixture filters
%particle.mixture_method = 'particles' ;
% Size of the different mixtures
particle.mixture_state_variables = 5 ;
particle.mixture_structural_shocks = 1 ;
particle.mixture_measurement_shocks = 1 ;
% Online approach
particle.liu_west_delta = 0.99 ;
particle.liu_west_max_resampling_tries = 5000;
% Options for setting the weights in conditional particle filters.
particle.cpf_weights_method.amisanotristani = true;
particle.cpf_weights_method.murrayjonesparslow = false;
particle.particle_filter_options ='';
% Copy particle structure in options_ global structure
options_.particle = particle;
options_.rwgmh.init_scale = 1e-4 ;
options_.rwgmh.scale_chain = 1 ;
options_.rwgmh.scale_shock = 1e-5 ;

% TeX output
options_.TeX = false;

% Exel
options_.xls_sheet = 1; % Octave does not support the empty string, rather use first sheet
options_.xls_range = '';

% Prior draws
options_.prior_draws = 10000;

% Prior posterior function sampling draws
options_.sampling_draws = 500;

options_.forecast = 0;
options_.forecasts.conf_sig = 0.9;
options_.conditional_forecast.conf_sig = 0.9;

% Model
options_.linear = false;

% Deterministic simulation
options_.stack_solve_algo = 0;
options_.markowitz = 0.5;
options_.minimal_solving_periods = 1;
options_.endogenous_terminal_period = false;
options_.no_homotopy = false;
options_.simul.endval_steady = false;

options_.simul.homotopy_max_completion_share = 1;
options_.simul.homotopy_min_step_size = 1e-3;
options_.simul.homotopy_step_size_increase_success_count = 3;
options_.simul.homotopy_initial_step_size = 1;
options_.simul.homotopy_linearization_fallback = false;
options_.simul.homotopy_marginal_linearization_fallback = 0; % Size of the step used for the marginal linearization; 0 means disabled

% Options used by perfect_foresight_* commands when they compute the steady
% state corresponding to a terminal condition
options_.simul.steady_solve_algo = options_.solve_algo;
options_.simul.steady_maxit = options_.steady.maxit;
options_.simul.steady_tolf = options_.solve_tolf;
options_.simul.steady_tolx = options_.solve_tolx;
options_.simul.steady_markowitz = options_.markowitz;

% Perfect foresight with expectation errors
options_.pfwee.constant_simulation_length = false;

% Solution
options_.order = 2;
options_.pruning = false;
options_.replic = 50;
options_.simul_replic = 1;
options_.drop = 100;
options_.aim_solver = false; % i.e. by default do not use G.Anderson's AIM solver, use mjdgges instead
options_.k_order_solver = false; % by default do not use k_order_perturbation but mjdgges
options_.partial_information = false;
options_.ACES_solver = false;
options_.conditional_variance_decomposition = [];

% Ramsey policy
options_.ramsey_policy = false;
options_.instruments = {};
options_.timeless = 0;
options_.ramsey.maxit = 500;
options_.ramsey.periods = 10000;
options_.ramsey.drop = 1000;

% estimation
options_.initial_period = NaN; %dates(1,1);
options_.no_init_estimation_check_first_obs=false;
options_.dataset.file = [];
options_.dataset.series = [];
options_.dataset.firstobs = dates();
options_.dataset.lastobs = dates();
options_.dataset.nobs = NaN;
options_.dataset.xls_sheet = [];
options_.dataset.xls_range = [];
options_.Harvey_scale_factor = 10;
options_.heteroskedastic_filter = false;
options_.MaxNumberOfBytes = 1e8;
options_.MaximumNumberOfMegaBytes = 111;
options_.analytic_derivation = 0; % Not a boolean, can also take values -1 or 2
options_.analytic_derivation_mode = 0;
options_.bayesian_irf = false;
options_.bayesian_th_moments = 0;
options_.diffuse_filter = false;
options_.filter_step_ahead = [];
options_.filtered_vars = false;
options_.smoothed_state_uncertainty = false;
options_.first_obs = NaN;
options_.nobs = NaN;
options_.kalman_algo = 0;
options_.fast_kalman_filter = false;
options_.kalman_tol = 1e-10;
options_.kalman.keep_kalman_algo_if_singularity_is_detected = false;
options_.diffuse_kalman_tol = 1e-6;
options_.use_univariate_filters_if_singularity_is_detected = 1;
options_.riccati_tol = 1e-6;
options_.lik_algo = 1;
options_.lik_init = 1;
options_.load_mh_file = false;
options_.load_results_after_load_mh = false;
options_.logdata = false;
options_.loglinear = false;
options_.linear_approximation = false;
options_.logged_steady_state = 0;
options_.mh_conf_sig = 0.90;
options_.prior_interval = 0.90;
options_.mh_drop = 0.5;
options_.mh_jscale = [];
options_.mh_tune_jscale.status = false;
options_.mh_tune_jscale.guess = [];
options_.mh_tune_jscale.target = .33;
options_.mh_tune_jscale.maxiter = 200000;
options_.mh_tune_jscale.rho = .7;
options_.mh_tune_jscale.stepsize = 1000;
options_.mh_tune_jscale.c1 = .02;
options_.mh_tune_jscale.c2 = .05;
options_.mh_tune_jscale.c3 = 4;
options_.mh_init_scale_factor = 2;
options_.mh_initialize_from_previous_mcmc.status = false;
options_.mh_initialize_from_previous_mcmc.directory = '';
options_.mh_initialize_from_previous_mcmc.record = '';
options_.mh_initialize_from_previous_mcmc.prior = '';
options_.mh_nblck = 2;
options_.mh_recover = false;
options_.mh_replic = 20000;
options_.recursive_estimation_restart = 0;
options_.MCMC_jumping_covariance='hessian';
options_.use_calibration_initialization = 0;
options_.endo_vars_for_moment_computations_in_estimation=[];
% occbin options
options_.occbin.likelihood.status=false;
options_.occbin.smoother.status=false;

% Run optimizer silently
options_.silent_optimizer = false;

% Prior restrictions
options_.prior_restrictions.status = 0;
options_.prior_restrictions.routine = [];

options_.mode_compute = 4;
options_.mode_file = '';
options_.moments_varendo = false;
options_.nk = 1;
options_.noconstant = false;
options_.nodiagnostic = false;
options_.mh_posterior_mode_estimation = false;
options_.prefilter = 0;
options_.presample = 0;
options_.prior_trunc = 1e-10;
options_.smoother = false;
options_.smoother_redux = false;
options_.posterior_max_subsample_draws = 1200;
options_.sub_draws = [];
options_.ME_plot_tol=1e-6;
options_.use_mh_covariance_matrix = false;
options_.gradient_method = 2; %used by csminwel and newrat
options_.gradient_epsilon = 1e-6; %used by csminwel and newrat
options_.posterior_sampler_options.sampling_opt = []; %extended set of options for individual posterior samplers
                                                      % Random Walk Metropolis-Hastings
options_.posterior_sampler_options.posterior_sampling_method = 'random_walk_metropolis_hastings';
options_.posterior_sampler_options.rwmh.proposal_distribution = 'rand_multivariate_normal';
options_.posterior_sampler_options.rwmh.student_degrees_of_freedom = 3;
options_.posterior_sampler_options.rwmh.use_mh_covariance_matrix=0;
options_.posterior_sampler_options.rwmh.save_tmp_file=0;
% Tailored Random Block Metropolis-Hastings
options_.posterior_sampler_options.tarb.proposal_distribution = 'rand_multivariate_normal';
options_.posterior_sampler_options.tarb.student_degrees_of_freedom = 3;
options_.posterior_sampler_options.tarb.mode_compute=4;
options_.posterior_sampler_options.tarb.new_block_probability=0.25; %probability that next parameter belongs to new block
options_.posterior_sampler_options.tarb.optim_opt=''; %probability that next parameter belongs to new block
options_.posterior_sampler_options.tarb.save_tmp_file=1;
% Slice
options_.posterior_sampler_options.slice.proposal_distribution = '';
options_.posterior_sampler_options.slice.rotated=0;
options_.posterior_sampler_options.slice.slice_initialize_with_mode=0;
options_.posterior_sampler_options.slice.use_mh_covariance_matrix=0;
options_.posterior_sampler_options.slice.WR=[];
options_.posterior_sampler_options.slice.mode_files=[];
options_.posterior_sampler_options.slice.mode=[];
options_.posterior_sampler_options.slice.initial_step_size=0.8;
options_.posterior_sampler_options.slice.save_tmp_file=1;
% Independent Metropolis-Hastings
options_.posterior_sampler_options.imh.proposal_distribution = 'rand_multivariate_normal';
options_.posterior_sampler_options.imh.use_mh_covariance_matrix=0;
options_.posterior_sampler_options.imh.save_tmp_file=0;
% Herbst and Schorfeide SMC Sampler
%options_.posterior_sampler = 'Herbst_Schorfheide' ;
options_.posterior_sampler_options.HSsmc.nphi= 25 ;
options_.posterior_sampler_options.HSsmc.lambda = 2 ;
options_.posterior_sampler_options.HSsmc.nparticles = 20000 ;
options_.posterior_sampler_options.HSsmc.c = 0.5 ;
options_.posterior_sampler_options.HSsmc.acpt = 0.25 ;
options_.posterior_sampler_options.HSsmc.trgt = 0.25 ;
options_.posterior_sampler_options.HSsmc.option_mutation = 1 ;
options_.posterior_sampler_options.HSsmc.alp = .9 ;
% DSMH: Dynamic Striated Metropolis-Hastings algorithm
%options_.posterior_sampler = 'DSMH' ;
options_.posterior_sampler_options.dsmh.H = 25 ;
options_.posterior_sampler_options.dsmh.N = 20 ;
options_.posterior_sampler_options.dsmh.G = 10 ;
options_.posterior_sampler_options.dsmh.K = 50 ;
options_.posterior_sampler_options.dsmh.lambda1 = 0.1 ;
options_.posterior_sampler_options.dsmh.nparticles = 20000 ;
options_.posterior_sampler_options.dsmh.alpha0 = 0.2 ;
options_.posterior_sampler_options.dsmh.alpha1 = 0.3 ; 
options_.posterior_sampler_options.dsmh.tau = 10 ;

options_.trace_plot_ma = 200;
options_.mh_autocorrelation_function_size = 30;
options_.plot_priors = 1;
options_.cova_compute = 1;
options_.parallel = 0;
options_.parallel_info.isHybridMatlabOctave = false;
options_.parallel_info.leaveSlaveOpen = 0;
options_.parallel_info.RemoteTmpFolder = '';
options_.parallel_info.use_psexec = true;
options_.number_of_grid_points_for_kde = 2^9;
quarter = 1;
years = [1 2 3 4 5 10 20 30 40 50];
options_.conditional_variance_decomposition_dates = zeros(1,length(years));
for i=1:length(years)
    options_.conditional_variance_decomposition_dates(i) = ...
        (years(i)-1)*4+quarter;
end
options_.filter_covariance = false;
options_.filter_decomposition = false;
options_.selected_variables_only = false;
options_.contemporaneous_correlation = false;
options_.initialize_estimated_parameters_with_the_prior_mode = 0;
options_.estimation.moments_posterior_density.indicator = true;
options_.estimation.moments_posterior_density.gridpoints = 2^9;
options_.estimation.moments_posterior_density.bandwidth = 0; % Rule of thumb optimal bandwidth parameter.
options_.estimation.moments_posterior_density.kernel_function = 'gaussian'; % Gaussian kernel for Fast Fourrier Transform approximaton.
                                                                            % Misc
                                                                            % options_.conf_sig = 0.6;

% homotopy for steady state
options_.homotopy_mode = 0;
options_.homotopy_steps = 1;
options_.homotopy_force_continue = false;

% numerical hessian
hessian.use_penalized_objective = false;

% Robust prediction error covariance (kalman filter)
options_.rescale_prediction_error_covariance = false;

options_.hessian = hessian;

%csminwel optimization routine
csminwel.tolerance.f=1e-7;
csminwel.maxiter=1000;
csminwel.verbosity=1;
csminwel.Save_files=false;

options_.csminwel=csminwel;

%newrat optimization routine
newrat.hess=1; % dynare numerical hessian
newrat.tolerance.f=1e-5;
newrat.tolerance.f_analytic=1e-7;
newrat.maxiter=1000;
newrat.verbosity=1;
newrat.Save_files=0;

options_.newrat=newrat;

% Simplex optimization routine (variation on Nelder Mead algorithm).
simplex.tolerance.x = 1e-4;
simplex.tolerance.f = 1e-4;
simplex.maxiter = 10000;
simplex.maxfcallfactor = 500;
simplex.maxfcall = [];
simplex.verbosity = 2;
simplex.delta_factor=0.05;
options_.simplex = simplex;

% CMAES optimization routine.
cmaes.SaveVariables='off';
cmaes.DispFinal='on';
cmaes.WarnOnEqualFunctionValues='no';
cmaes.DispModulo='10';
cmaes.LogModulo='0';
cmaes.LogTime='0';
cmaes.TolFun = 1e-7;
cmaes.TolX = 1e-7;
cmaes.Resume = 0;
options_.cmaes = cmaes;

% simpsa optimization routine.
simpsa.TOLFUN = 1e-4;
simpsa.TOLX = 1e-4;
simpsa.TEMP_END = .1;
simpsa.COOL_RATE = 10;
simpsa.INITIAL_ACCEPTANCE_RATIO = .95;
simpsa.MIN_COOLING_FACTOR = .9;
simpsa.MAX_ITER_TEMP_FIRST = 50;
simpsa.MAX_ITER_TEMP_LAST = 2000;
simpsa.MAX_ITER_TEMP = 10;
simpsa.MAX_ITER_TOTAL = 5000;
simpsa.MAX_TIME = 2500;
simpsa.MAX_FUN_EVALS = 20000;
simpsa.DISPLAY = 'iter';
options_.simpsa = simpsa;

%solveopt optimizer
solveopt.minimizer_indicator=-1; %use minimizer
solveopt.TolX=1e-6; %accuracy of argument
solveopt.TolFun=1e-6; %accuracy of function
solveopt.MaxIter=15000;
solveopt.verbosity=1;
solveopt.TolXConstraint=1.e-8;
solveopt.SpaceDilation=2.5;
solveopt.LBGradientStep=1.e-11;
options_.solveopt=solveopt;

%simulated annealing
options_.saopt.neps=10;
options_.saopt.maximizer_indicator=0;
options_.saopt.rt=0.10;
options_.saopt.MaxIter=100000;
options_.saopt.verbosity=1;
options_.saopt.TolFun=1.0e-8;
options_.saopt.initial_temperature=15;
options_.saopt.ns=10;
options_.saopt.nt=10;
options_.saopt.step_length_c=0.1;
options_.saopt.initial_step_length=1;

% particleswarm (global optimization toolbox needed)
particleswarm.Display = 'iter';
particleswarm.DisplayInterval = 1;
particleswarm.FunctionTolerance = 1e-6;
particleswarm.FunValCheck = 'on';
particleswarm.HybridFcn = [];
particleswarm.InertiaRange = [0.1, 1.1];
particleswarm.MaxIterations = 100000;
particleswarm.MaxStallIterations = 20;
particleswarm.MaxStallTime = Inf;
particleswarm.MaxTime = Inf;
particleswarm.MinNeighborsFraction = .25;
particleswarm.ObjectiveLimit = -Inf;
particleswarm.UseParallel = false;
particleswarm.UseVectorized = false;
options_.particleswarm = particleswarm;

% prior analysis
options_.prior_mc = 20000;
options_.prior_analysis_endo_var_list = {};

% did model undergo block decomposition + minimum feedback set computation ?
options_.block = false;

% model evaluated using a compiled MEX
options_.use_dll = false;

% model evaluated using bytecode.dll
options_.bytecode = false;

% if true, use a fixed point method to solve Sylvester equation (for large scale models)
options_.sylvester_fp = false;

% convergence criteria to solve iteratively a sylvester equations
options_.sylvester_fixed_point_tol = 1e-12;

% if true, use a fixed point method to solve Lyapunov equation (for large scale models)
options_.lyapunov_fp = false;
% if true, use a doubling algorithm to solve Lyapunov equation (for large scale models)
options_.lyapunov_db = false;
% if true, use a square root solver to solve Lyapunov equation (for large scale models)
options_.lyapunov_srs = false;

% convergence criterion for iteratives methods to solve lyapunov equations
options_.lyapunov_fixed_point_tol = 1e-10;
options_.lyapunov_doubling_tol = 1e-16;

% if true, use a cycle reduction method to compute the decision rule (for large scale models)
options_.dr_cycle_reduction = false;

% convergence criterion for iteratives methods to solve the decision rule
options_.dr_cycle_reduction_tol = 1e-7;

% if true, use a logarithmic reduction method to compute the decision rule (for large scale models)
options_.dr_logarithmic_reduction = false;

% convergence criterion for iteratives methods to solve the decision rule
options_.dr_logarithmic_reduction_tol = 1e-12;

% convergence criterion for iteratives methods to solve the decision rule
options_.dr_logarithmic_reduction_maxiter = 100;

% dates for historical time series
options_.initial_date = dates();

% discretionary policy
options_.discretionary_policy = false;
options_.discretionary_tol = 1e-7;

% Shock decomposition
options_.parameter_set = [];
options_.use_shock_groups = '';
options_.shock_decomp.colormap = '';
options_.shock_decomp.init_state = 0;
options_.shock_decomp.with_epilogue = false;

% Shock decomposition realtime
options_.shock_decomp.forecast = 0;
options_.shock_decomp.presample = NaN;
options_.shock_decomp.save_realtime = 0; % saves memory
options_ = set_default_plot_shock_decomposition_options(options_);

% Nonlinearfilters
options_.nonlinear_filter = [];

% SBVAR & MS SBVAR initializations:
% SBVAR
options_.ms.vlistlog = [];
options_.ms.restriction_fname = 0;
options_.ms.cross_restrictions = false;
options_.ms.contemp_reduced_form = false;
options_.ms.real_pseudo_forecast = 0;
options_.ms.dummy_obs = 0;
options_.ms.ncsk = 0;
options_.ms.indxgforhat = 1;
options_.ms.indxgimfhat = 1;
options_.ms.indxestima = 1;
options_.ms.indxgdls = 1;
options_.ms.cms =0;
options_.ms.ncms = 0;
options_.ms.eq_cms = 1;
options_.ms.banact = 1;
options_.ms.log_var = [];
options_.ms.Qi = [];
options_.ms.Ri = [];
options_.ms.lower_cholesky = 0;
options_.ms.upper_cholesky = 0;
options_.ms.constants_exclusion = 0;

% MS SBVAR (and some SBVAR)
options_ = initialize_ms_sbvar_options(M_, options_);

% saved graph formats
options_.graph_save_formats.eps = 1;
options_.graph_save_formats.pdf = 0;
options_.graph_save_formats.fig = 0;

% endogenous prior
options_.endogenous_prior = false;
options_.endogenous_prior_restrictions.irf={};
options_.endogenous_prior_restrictions.moment={};

% OSR Optimal Simple Rules
options_.osr.opt_algo=4;

% use GPU
options_.gpu = false;

%Geweke convergence diagnostics
options_.convergence.geweke.taper_steps=[4 8 15];
options_.convergence.geweke.geweke_interval=[0.2 0.5];
%Raftery/Lewis convergence diagnostics;
options_.convergence.rafterylewis.indicator=false;
options_.convergence.rafterylewis.qrs=[0.025 0.005 0.95];

%tolerance for Modified Harmonic Mean estimator
options_.marginal_data_density.harmonic_mean.tolerance = 0.01;

% Options for lmmcp solver
options_.lmmcp.status = false;

% Options for lcppath solver
options_.lcppath.A = [];
options_.lcppath.b = [];
options_.lcppath.t = [];
options_.lcppath.mu0 = [];

% Options for mcppath solver
options_.mcppath.A = [];
options_.mcppath.b = [];
options_.mcppath.t = [];
options_.mcppath.mu0 = [];

%Figure options
options_.figures.textwidth=0.8;

options_.varobs_id=[]; %initialize field

options_.pac.estimation.ols.share_of_optimizing_agents.lb = 0.0;
options_.pac.estimation.ols.share_of_optimizing_agents.ub = 1.0;

options_.conditional_likelihood.status = false;
options_.conditional_likelihood.order = 1;
