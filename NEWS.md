Announcement for Dynare 5.2 (on 2022-07-27)
===========================================

We are pleased to announce the release of Dynare 5.2.

This maintenance release fixes various bugs.

The Windows, macOS and source packages are already available for download at
[the Dynare website](https://www.dynare.org/download/).

All users are strongly encouraged to upgrade.

This release is compatible with MATLAB versions ranging from 8.3 (R2014a) to
9.12 (R2022a), and with GNU Octave version 6.4.0 (under Windows).

Note for macOS users with an Apple Silicon processor, and who are also MATLAB
users: the official MATLAB version for use on those processors is still the
Intel version (running through Rosetta 2), so the official Dynare package for
download on our website is built for Intel only. However, since Mathworks has
released a beta version of MATLAB for Apple Silicon, we created a beta package
of Dynare that you can try with it. See this forum thread for more details:
https://forum.dynare.org/t/testers-wanted-release-of-dynare-5-x-beta-for-apple-silicon-m1-m2-chips/20499

Here is a list of the problems identified in version 5.1 and that have been
fixed in version 5.2:

* Problems with the `steady_state` operator:
  + if a `steady_state` operator contained an algebraic expression appearing
    multiple times in the model and sufficiently complex to trigger the
    creation of a temporary term, then the result of the operator would be
    wrong (the operator was essentially ignored)
  + if a `steady_state` operator contained a call to an external function, then
    the result of the operator would be wrong (the operator was essentially
    ignored). A proper fix to this problem would require substantial
    architectural changes, so for now it is forbidden to use an external
    function inside a `steady_state` operator
* Pruning in particle filtering at order 2 was not using the exact same formula
  as the original Kim et al. (2008) paper. A second-order term entered the
  cross-product between states and shocks, where it should have been a
  first-order term. This however would not lead to explosive trajectories in
  practice
* The `simul_replic` option of the `stoch_simul` command would not store the
  binary file in the `Output` folder
* Problems with Ramsey policy (`ramsey_model`/`ramsey_policy` commands):
  + steady state files would not work when auxiliary variables included
    Lagrange multipliers
  + for linear competitive equilibrium laws of motion, welfare evaluated at
    higher order was erroneously equated to steady state welfare
* The `discretionary_policy` command would not always correctly infer the
  number of instruments and equations, leading to spurious error messages
* Perfect foresight simulations of purely forward or backward models would
  crash if complex numbers were encountered
* When using both `block` and `bytecode` options of the `model` block, if the
  model was such that a sufficiently complex algebraic expression appeared both
  in the residuals and in the derivatives, leading to the creation of a
  temporary term, then the results could be incorrect under some circumstances
* When using the `bytecode` option of the `model` block, leads of more than
  +127 or lags of less than -128 were not correctly handled
* Problems with the solver under occasionally binding constraints
  (`occbin_solver` command):
  + when solving the baseline regime, it would not properly handle errors like
    Blanchard-Kahn violations
  + the piecewise linear Kalman filter (PKF) would crash if the model solution
    could not be computed for a parameter draw
  + the `oo_.FilteredVariablesKStepAhead` and `oo_.UpdatedVariables`
    MATLAB/Octave variables would contain the steady state twice
  + the inversion filter would crash if the `filter_step_ahead` or
    `state_uncertainty` options were requested
  + the PKF would crash if `filter_step_ahead=1` was specified
  + the PKF would crash if the `state_uncertainty` option was specified
    together with the `smoother_redux` option
  + the last regime before the system is back to normal times in the
    two-constraints case could be wrongly set, possibly leading to wrong
    simulations, lack of convergence or crashes
* Problems with identification (`identification` command):
  + with `prior_mc>1` specified, it would incorrectly display the share of rank
    deficient Jacobians
  + it would crash during plotting or displaying identification strength when
    the necessary identification criteria based on moments could not be
    computed
* The `plot_shock_decomposition` command would crash if invalid field names
  were encountered
* The `shock_decomposition` command would not pass specified initial dates to
  generated plots
* Various pathological cases encountered in steady state finding could lead to
  a crash
* The `solve_algo=0` option of the `steady` command would not honor `tolx`
* In the `dynare_sensitivity` command, stability mapping would not correctly
  honor the prior bounds


Announcement for Dynare 5.1 (on 2022-04-06)
===========================================

We are pleased to announce the release of Dynare 5.1.

This maintenance release fixes various bugs.

The Windows, macOS and source packages are already available for download at
[the Dynare website](https://www.dynare.org/download/).

All users are strongly encouraged to upgrade.

This release is compatible with MATLAB versions ranging from 8.3 (R2014a) to
9.12 (R2022a), and with GNU Octave version 6.4.0 (under Windows).

Here is a list of the problems identified in version 5.0 and that have been
fixed in version 5.1:

* Various problems with perfect foresight simulations in combination with
  `block` and/or `bytecode` options of the `model` block:
  + Simulation with `bytecode` and `stack_solve_algo=4` could give incorrect
    results if the model has a linear block of type “Solve two boundaries
    simple/complete”
  + Simulation with `bytecode` and `stack_solve_algo=1` could fail to converge
  + Simulation with `block` (but without `bytecode`) and `stack_solve_algo=1`
    gave wrong results in the last simulation period if the model has a block
    of type “Solve two boundaries simple/complete”
  + Simulation with `bytecode` and `block` would give incorrect results if the
    model has a linear block of type “Solve forward simple/complete”
  + Simulation with `block` (but without `bytecode`) would crash or give
    incorrect results if the model has a block of type “Solve forward/backward
    simple/complete”
  + Simulation with `bytecode`, `block` and `stack_solve_algo={0,1,4}` would
    crash or give incorrect results if the model has a block of type “Solve
    forward/backward complete”
  + Simulation with `block` (but without `bytecode`) gave incorrect results if
    the model has a block of type “Solve backward simple/complete”
  + Simulation with `block` (with or without `bytecode`) could give incorrect
    results if the model has a nonlinear block of type “Solve forward/backward
    simple/complete”
  + Simulation with `bytecode`, `block` and `stack_solve_algo=4` could give
    incorrect results if the model has a block of type “Solve backward/forward
    simple/complete” that follows a block of type “Solve two boundaries” (in
    the sense of the dependency graph)
  + The convergence criterion in simulations with `block` (but without
    `bytecode`) was incorrect: the value of the `tolf` option from the `steady`
    command was used instead of the value of `tolf` option from the
    `perfect_foresight_solver` command
* Various problems with steady state computation in combination with `block`
  and/or `bytecode` options of the `model` block:
  + Steady state computation with `bytecode` and `block` could fail if some
    equations are marked `[static]`
  + Steady state computation with `bytecode`, `block` and `solve_algo` ⩽ 4 or ⩾
    9 could fail
  + Steady state computation with `bytecode`, `block` and `solve_algo=6` would
    crash or give incorrect results if the model has a block of type “Solve
    forward/backward complete”
* The `check` command would crash or give incorrect results when using the
  `block` option of the `model` block and if the model has a block of type
  “Solve backward complete”
* The `static` and `incidence` options of the `model_info` command did not work
  as documented in the reference manual
* Various problems with the `method_of_moments` command:
  + It would crash if no `matched_moments` block is present
  + It would always load the full range of the first Excel sheet instead of the
    `xls_range` of the specified `xls_sheet`
  + SMM would crash if a parameter draw triggers an error during
    `additional_optimizer_steps = 13`
  + The `debug` option could not be passed to the command
* In the `estimation` command, the `scale_file` field of the
  `posterior_sampler_options` option did not correctly load the scale
* The `moments_varendo` option of the `estimation` command could crash for
  large models
* The `resid` command would not show `name` tags when used in conjunction with
  the `ramsey_model` command
* Simulations with the `occbin_solver` command would not work if there is only
  a surprise shock in the first period
* The Liu & West auxiliary particle filter could enter infinite loops


Announcement for Dynare 5.0 (on 2022-01-07)
===========================================

We are pleased to announce the release of Dynare 5.0.

This major release adds new features and fixes various bugs.

The Windows, macOS and source packages are already available for download at
[the Dynare website](https://www.dynare.org/download/).

All users are strongly encouraged to upgrade.

This release is compatible with MATLAB versions ranging from 8.3 (R2014a) to
9.11 (R2021b), and with GNU Octave version 6.4.0 (under Windows).

The new tools for semi-structural models and the improvements on the nonlinear
solvers were funded by the European Central Bank. Special thanks to Nikola
Bokan (ECB) for his contributions and numerous bug reports and fixes.

Major user-visible changes
--------------------------

 - New routines for simulating semi-structural (backward) models where
   some equations incorporate expectations based on future values of a VAR or
   trend component model. See the `var_model`, `trend_component_model` and
   `var_expectation_model` commands, and the `var_expectation` operator.

 - New routines for simulating semi-structural models where some equations are
   specified using the polynomial adjustment costs (PAC) approach, as in the
   FRB/US model (see Brayton et al., 2014 and Brayton et al., 2000) and the
   ECB-BASE model (see Angelini et al., 2019). The forward-looking terms of the
   PAC equations can be computed either using a satellite VAR model, or using
   full model-consistent expectations. See the `pac_model` command and the
   `pac_expectation` operator.

 - New Method of Moments toolbox that provides functionality to estimate
   parameters by (i) Generalized Method of Moments (GMM) up to 3rd-order pruned
   perturbation approximation or (ii) Simulated Method of Moments (SMM) up to
   any perturbation approximation order. The toolbox is inspired by replication
   codes accompanying Andreasen et al. (2018), Born and Pfeifer (2014), and
   Mutschler (2018). It is accessible via the new `method_of_moments` command
   and the new `matched_moments` block. Moreover, by default, a new non-linear
   least squares optimizer based on `lsqnonlin` is used for minimizing the
   method of moments objective function (available under `mode_compute=13`).
   GMM can further benefit from using gradient-based optimizers (using
   `analytic_standard_errors` option and/or passing `'Jacobian','on'` to the
   optimization options) as the Jacobian of the moment conditions can be
   computed analytically.

 - Implementation of the Occbin algorithm by Guerrieri and Iacoviello (2015),
   together with the inversion filter of Cuba-Borda, Guerrieri, Iacoviello, and
   Zhong (2019) and the piecewise Kalman filter of Giovannini, Pfeiffer, and
   Ratto (2021). It is available via the new block `occbin_constraints` and the
   new commands `occbin_setup`, `occbin_solver`, `occbin_graph`, and
   `occbin_write_regimes`.

 - Stochastic simulations

    - `stoch_simul` now supports theoretical moments at `order=3` with
      `pruning`.

    - `stoch_simul` now reports second moments based on the pruned state space
      if the `pruning` option is set (in previous Dynare releases it would
      report a second-order accurate result based on the linear solution).

 - Estimation

    - Performance optimization to pruned state space systems and Lyapunov
      solvers.

    - New option `mh_posterior_mode_estimation` to `estimation` to perform
      mode-finding by running the MCMC.

    - New heteroskedastic filter and smoother, where shock standard errors may
      *unexpectedly* change in every period. Triggered by the
      `heteroskedastic_filter` option of the `estimation` command, and
      configured via the `heteroskedastic_shocks` block.

    - New option `mh_tune_guess` for setting the initial value for
      `mh_tune_jscale`.

    - New option `smoother_redux` to `estimation` and `calib_smoother` to
      trigger computing the Kalman smoother on a restricted state space instead
      of the full one.

    - New block `filter_initial_state` for setting the initial condition of the
      Kalman filter/smoother.

    - New option `mh_initialize_from_previous_mcmc` to the `estimation` command
      that allows to pick initial values for a new MCMC from a previous one.

    - The `xls_sheet` option of the `estimation` command now takes a quoted
      string as value. The former unquoted syntax is still accepted, but no
      longer recommended.

    - New option `particle_filter_options` to set various particle filter options.

 - Perfect foresight and extended path

    - New specialized algorithm in `perfect_foresight_solver` to deal with
      purely static problems.

    - The `debug` option of `perfect_foresight_solver` provides debugging
      information if the Jacobian is singular.

    - In deterministic models (perfect foresight or extended path), exogenous
      variables with lead/lags are now replaced by auxiliary variables. This
      brings those models in line with the transformation done on stochastic
      models. However, note that the transformation is still not exactly the same
      between the two classes of models, because there is no need to take into
      account the Jensen inequality for the latter. In deterministic models,
      there is a one-to-one mapping between exogenous with lead/lags and
      auxiliaries, while in stochastic models, an auxiliary endogenous may
      correspond to a more complex nonlinear expression.

 - Optimal policy

    - Several improvements to `evaluate_planner_objective`:

       - it now applies a consistent approximation order when doing the
         computation;
       - in addition to the conditional welfare, it now also provides the
         unconditional welfare;
       - in a stochastic context, it now works with higher order approximation
         (only the conditional welfare is available for order ⩾ 3);
       - it now also works in a perfect foresight context.

    - `discretionary_policy` is now able to solve nonlinear models (it will
      then use their first-order approximation, and the analytical steady state
      must be provided).

 - Identification

    - New option `schur_vec_tol` to the `identification` command, for setting
      the tolerance level used to find nonstationary variables in the Schur
      decomposition of the transition matrix.

    - The `identification` command now supports optimal policy.

 - Shock decomposition

    - The `fast_realtime` option of the `realtime_shock_decomposition` command
      now accepts a vector of integers, which runs the smoother for all the
      specified data vintages.

 - Macro processor

    - Macroprocessor variables can be defined without a value (they are
      assigned integer 1).

 - LaTeX and JSON outputs

    - New `nocommutativity` option to the `dynare` command. This option tells
      the preprocessor not to use the commutativity of addition and
      multiplication when looking for common subexpressions. As a consequence,
      when using this option, equations in various outputs (LaTeX, JSON…) will
      appear as the user entered them (without terms or factors swapped). Note
      that using this option may have a performance impact on the preprocessing
      stage, though it is likely to be small.

    - Model-local variables are now substituted out as part of the various
      model transformations. This means that they will no longer appear in
      LaTeX or in JSON files (for the latter, they are still visible with
      `json=parse` or `json=check`).

 - Compilation of the model (`use_dll` option)

    - Block decomposition (option `block` of `model`) can now be used in
      conjunction with the `use_dll` option.

    - The `use_dll` option can now directly be given to the `dynare` command.

 - dseries classes

    - Routines for converting between time series frequencies (e.g. daily to
      monthly) have been added.

    - dseries now supports bi-annual and daily frequency data.

    - dseries can now import data from [DBnomics](https://db.nomics.world), via
      the [mdbnomics](https://git.dynare.org/dbnomics/mdbnomics) plugin. Note
      that this does not yet work under Octave. For the time being, the
      DBnomics plugin must be installed separately.

 - Misc improvements

    - The `histval_file` and `initval_file` commands have been made more
      flexible and now have functionalities similar to the `datafile` option of
      the `estimation` command.

    - When using the `loglinear` option, the output from Dynare now clearly
      shows that the results reported concern the log of the original variable.

    - Options `block` and `bytecode` of `model` can now be used in conjunction
      with model-local variables (variables declared with a pound-sign `#`).

    - The `model_info` command now prints the typology of endogenous variables
      for non-block decomposed models.

    - The total computing time of a run (in seconds) is now saved to `oo_.time`.

    - New `notime` option to the `dynare` command, to disable the printing and
      the saving of the total computing time.

    - New `parallel_use_psexec` command-line Windows-specific option for
      parallel local clusters: when `true` (the default), use `psexec` to spawn
      processes; when `false`, use `start`.

    - When compiling from source, it is no longer necessary to pass the
      `MATLAB_VERSION` version to the configure script; the version is now
      automatically detected.

Incompatible changes
--------------------

 - Dynare will now generally save its output in the `MODFILENAME/Output` folder
   (or the `DIRNAME/Output` folder if the `dirname` option was specified)
   instead of the main directory. Most importantly, this concerns the
   `_results.mat` and the `_mode.mat` files.

 - The structure of the `oo_.planner_objective` field has been changed, in
   relation to the improvements to `evaluate_planner_objective`.

 - The preprocessor binary has been renamed to `dynare-preprocessor`, and is
   now located in a dedicated `preprocessor` subdirectory.

 - The `dynare` command no longer accepts `output=dynamic` and `output=first`
   (these options actually had no effect).

 - The minimal required MATLAB version is now R2014a (8.3).

 - The 32-bit support has been dropped for Windows.

Bugs that were present in 4.6.4 and that have been fixed in 5.0
---------------------------------------------------------------

* Equations marked with `static`-tags were not detrended when a `deflator` was
  specified
* Parallel execution of `dsge_var` estimation was broken
* The preprocessor would incorrectly simplify forward-looking constant
  equations of the form `x(+1)=0` to imply `x=0`
* Under some circumstances, the use of the `model_local_variable` statement
  would lead to a crash of the preprocessor
* When using the `block`-option without `bytecode` the residuals of the static
  model were incorrectly displayed
* When using `k_order_solver`, the `simult_` function ignored requested
  approximation orders that differed from the one used to compute the decision
  rules
* Stochastic simulations of the `k_order_solver` without `pruning` iterated on
  the policy function with a zero shock vector for the first (non-endogenous)
  period
* `estimation` would ignore the mean of non-zero observables if the mean was 0
  for the initial parameter vector
* `mode_check` would crash if a parameter was estimated to be exactly 0
* `load_mh_file` would not be able to load the proposal density if the previous run
  was done in parallel
* `load_mh_file` would not work with MCMC runs from Dynare versions before
  4.6.2
* `ramsey_model` would not correctly work with `lmmcp`
* `ramsey_model` would crash if a non-scalar error code was encountered during
  steady state finding.
* Using undefined objects in the `planner_objective` function would yield an
  erroneous error message about the objective containing exogenous variables
* `model_diagnostics` did not correctly handle a previous `loglinear` option
* `solve_algo=3` (csolve) would ignore user-set `maxit` and `tolf` options
* The `planner_objective` values were not based on the correct initialization
  of auxiliary variables (if any were present)
* The `nostrict` command line option was not ignoring unused endogenous
  variables in `initval`, `endval`, and `histval`
* `prior_posterior_statistics_core` could crash for models with eigenvalues
  very close to 1
* The display of the equation numbers in `debug` mode related to issues in the
  Jacobian would not correctly take auxiliary equations into account
* The `resid` command was not correctly taking auxiliary and missing equations
  related to optimal policy (`ramsey_model`, `discretionary_policy`) into
  account
* `bytecode` would lock the `dynamic.bin` file upon encountering an exception,
  requiring a restart of MATLAB to be able to rerun the file
* Estimation with the `block` model option would crash when calling the block
  Kalman filter
* The `block` model option would crash if no `initval` statement was present
* Having a variable with the same name as the mod-file present in the base
  workspace would result in a crash
* `oo_.FilteredVariablesKStepAheadVariances` was wrongly computed in the Kalman
  smoother based on the previous period forecast error variance
* Forecasts after `estimation` would not work if there were lagged exogenous
  variables present
* Forecasts after `estimation` with MC would crash if measurement errors were
  present
* Smoother results would be infinity for auxiliary variables associated with
  lagged exogenous variables
* In rare cases, the posterior Kalman smoother could crash due to previously
  accepted draws violating the Blanchard-Kahn conditions when using an
  unrestricted state space
* `perfect_foresight_solver` would crash for purely static problems
* Monte Carlo sampling in `identification` would crash if the minimal state
  space for the Komunjer and Ng test could not be computed
* Monte Carlo sampling in `identification` would skip the computation of
  identification statistics for all subsequent parameter draws if an error was
  triggered by one draw
* The `--steps`-option of Dynare++ was broken
* `smoother2histval` would crash if variable names were too similar
* `smoother2histval` was not keeping track of whether previously stored results
  were generated with `loglinear`
* The `initval_file` option was not supporting Dynare’s translation of a model
  into a one lead/lag-model via auxiliary variables

References
----------

 - Andreasen et al. (2018): “The pruned state-space system for non-linear DSGE
   models: Theory and empirical applications,” Review of Economic Studies,
   85(1), 1–49

 - Angelini, Bokan, Christoffel, Ciccarelli and Zimic (2019): “Introducing
   ECB-BASE: The blueprint the new ECB semi-structural model for the euro area,”
   ECB Working Paper no. 2315

 - Born and Pfeifer (2014): “Policy risk and the business cycle,” Journal of
   Monetary Economics, 68, 68–85

 - Brayton, Davis and Tulip (2000): “Polynomial adjustment costs in FRB/US,”
   Unpublished manuscript

 - Brayton, Laubach, and Reifschneider (2014): “The FRB/US Model: A tool for
   macroeconomic policy analysis,” FEDS Notes. Washington: Board of Governors
   of the Federal Reserve System, https://doi.org/10.17016/2380-7172.0012

 - Cuba-Borda, Guerrieri, Iacoviello, and Zhong (2019): “Likelihood evaluation
   of models with occasionally binding constraints,” Journal of Applied
   Econometrics, 34(7), 1073–1085

 - Giovannini, Pfeiffer, and Ratto (2021): “Efficient and robust inference of
   models with occasionally binding constraints,” Working Paper 2021-03, Joint
   Research Centre, European Commission

 - Guerrieri and Iacoviello (2015): “OccBin: A toolkit for solving dynamic
   models with occasionally binding constraints easily,” Journal of Monetary
   Economics, 70, 22–38

 - Mutschler (2018): “Higher-order statistics for DSGE models,” Econometrics
   and Statistics, 6(C), 44–56


Announcement for Dynare 4.6.4 (on 2021-03-18)
=============================================

We are pleased to announce the release of Dynare 4.6.4.

This maintenance release fixes various bugs.

The Windows, macOS and source packages are already available for download at
[the Dynare website](https://www.dynare.org/download/).

All users are strongly encouraged to upgrade.

This release is compatible with MATLAB versions ranging from 7.9 (R2009b) to
9.10 (R2021a), and with GNU Octave version 6.2.0 (under Windows).

Here is a list of the problems identified in version 4.6.3 and that have been
fixed in version 4.6.4:

* Passing multiple shock values through a MATLAB/Octave vector in a `mshocks`
  block would not work
* The `mode_compute=12` option was broken
* The `use_mh_covariance_matrix` option was ignored
* The `load_mh_file` option together with `mh_replic=0` would not allow
  computing `moments_varendo` for a different list of variables
* The `forecast` option was ignored when using `mode_compute=0` without a
  mode-file to execute the smoother
* The `discretionary_policy` command would crash in the presence of news shocks
* The `ramsey_constraints` block would crash if the constraints contained
  defined `parameters`
* Identification would display a wrong error message if a unit root was present
  and `diffuse_filter` had been set
* Particle filter estimation would crash if the initial state covariance became
  singular for a draw
* Particle filter estimation would crash if `k_order_solver` option was
  specified with `options_.particle.pruning`
* The initial state covariance in particle filter estimation could be `NaN`
  when using `nonlinear_filter_initialization=2` despite
  `options_.particles.pruning=1`
* Output of `smoother` results when using particle filters would be based on
  `order=1`
* Output of `moments_varendo` results when using particle filters would be
  based on `order=1`
* When decreasing the `order` in `.mod` files, `oo_.dr` could contain stale
  results from higher orders
* Estimation results using the particle filter at `order=3` would be incorrect
  if the restricted state space differed from the unrestricted one
* The `mode_compute=102` option (SOLVEOPT) could return with `Inf` instead of
  the last feasible value
* Using `analytic_derivation` for Bayesian estimation would result in wrong
  results when the multivariate Kalman filter entered the steady state stage
* Using `analytic_derivation` for maximum likelihood estimation would result in
  a crash
* When using the Bayesian smoother with `filtered_vars`, the field for
  `Filtered_Variables_X_step_ahead` used the length of vector instead of the
  actual steps in `filter_step_ahead`
* `mode_compute=1,3` crashed when `analytic_derivation` was specified
* `mode_compute=1,3,102` did only allow for post-MATLAB 2016a option names
* The `cova_compute=0` option was not working with user-defined
  `MCMC_jumping_covariance`
* The `mode_compute=1` option was not working with `analytic_derivation`
* Not all commands were honouring the `M_.dname` folder when saving
* LaTeX output of the simulated variance decomposition for observables with
  measurement error could have a wrong variable label

Announcement for Dynare 4.6.3 (on 2020-11-23)
=============================================

We are pleased to announce the release of Dynare 4.6.3.

This maintenance release fixes various bugs.

The Windows, macOS and source packages are already available for download at
[the Dynare website](https://www.dynare.org/download/).

All users are strongly encouraged to upgrade.

This release is compatible with MATLAB versions ranging from 7.9 (R2009b) to
9.9 (R2020b), and with GNU Octave versions 5.2.0 (under Windows) and 4.4.1
(under macOS).

Here is a list of the problems identified in version 4.6.2 and that have been
fixed in version 4.6.3:

* Using an unknown symbol in `irf_shocks` option of `stoch_simul` would lead to
  a crash of the preprocessor
* `identification` would crash for purely forward-looking models
* The `endogenous_prior` option did not properly handle missing observations
* The auxiliary particle filter with pruning and resampling would crash
* Initialization of the state variance for particle filters was buggy
* An `@#else` clause after an `@#ifndef` was not correctly interpreted
* An `@#elseif` clause after an `@#ifdef` or an `@#ifndef` was not correctly
  interpreted
* Perfect foresight simulations of models with a single equation would crash
  when using either the `lmmcp` option or the `linear_approximation`
* Inequality constraints on endogenous variables (when using the `lmmcp`
  option) were not enforced on purely backward or purely forward models
* Perfect foresight simulations with `bytecode` and `block` options could crash
  if there was a purely forward variable whose value in all periods could be
  evaluated backward (typically a process of the form `y=a*y(+1)+e`)
* `extended_path` was broken with `bytecode`
* Under Windows, with Octave, the k-order perturbation and MS-SBVAR MEX files
  could not be loaded
* On Fedora (and possibly other GNU/Linux distributions), compilation from
  source would fail against Octave 5

Announcement for Dynare 4.6.2 (on 2020-09-07)
=============================================

We are pleased to announce the release of Dynare 4.6.2.

This maintenance release fixes various bugs.

The Windows, macOS and source packages are already available for download at
[the Dynare website](https://www.dynare.org/download/).

All users are strongly encouraged to upgrade.

This release is compatible with MATLAB versions ranging from 7.9 (R2009b) to
9.8 (R2020a), and with GNU Octave versions 5.2.0 (under Windows) and 4.4.1
(under macOS).

*Note for Windows users:* upon launching the Dynare installer, you may get a
warning emitted by Windows Defender SmartScreen, saying that this is an
unrecognized app and that it was prevented from starting. You can safely ignore
this warning, as long as you can verify on the next screen that CEPREMAP is the
editor of the software. This security warning is due to the fact that we had to
renew our code signing certificate (which had expired), and it takes some time
to rebuild our reputation as a software editor using the new certificate.

Here is a list of the problems identified in version 4.6.1 and that have been
fixed in version 4.6.2:

* Perfect foresight simulations of purely backward models could deliver an
  incorrect result if some exogenous variable appeared with a lag of 2 or more
  (and neither `block` nor `bytecode` option was used)
* Perfect foresight simulations of linear models could deliver an incorrect
  result if the following four conditions were met:
  + the model was actually declared as linear through the `linear` option
  + there was an exogenous variable with a lead or a lag
  + `stack_solve_algo` was equal to 0 (the default) or 7
  + neither `block` nor `bytecode` option was used
* In stochastic simulations, for variables that actually do not leave the
  steady state, reported simulated moments could be spurious (due to division
  by zero)
* Displayed variance decompositions would only take into account measurement
  errors if measurement errors were present for all observed variables
* The posterior variance decompositions with measurement errors computed with
  `moments_varendo` were incorrect
* `moments_varendo` would not update `oo_.PosteriorTheoreticalMoments` if it
  was already present, from *e.g.* an earlier run of `estimation`
* Identification would in some cases compute wrong Jacobian of moments
* Identification would display incorrect results if parameter dependence was
  implemented via a steady state file
* `generate_trace_plots` would crash when measurement errors were present
* `estimation` would crash for correlated measurement errors
* Parallel execution/testing could crash instead of aborting with a proper
  error message
* Under macOS, Dynare would incorrectly claim that it is compiled for Octave
  5.2.0 (it is actually compiled for Octave 4.4.1)
* Using external functions in a model local variable would crash the
  preprocessor
* Tolerance criteria for steady state computations were inconsistently set
* `stoch_simul` with its default `order=2` would crash with a message about
  `hessian_eq_zero` not existing if an explicit `order=1` was present somewhere
  else in the `.mod` file
* Model local variables were not written to the `modfile.json` JSON file
* Model local variables names would have two spurious underscores at their
  point of definition in the `dynamic.json` and `static.json` files (but only
  in the definition, not when they were used, which is inconsistent)
* The `solve_algo=9` option was not accessible. The `solve_algo=10` and
  `solve_algo=11` options were not accessible with `block` (without `bytecode`)
* Under certain circumstances, `extended_path` would crash when used in
  conjunction with the `block` option
* `extended_path` was not working with the `bytecode` option
* `shock_decomposition` was not accepting the options of `estimation` related
  to smoothing
* `conditional_forecast` would display a warning even if the simulation was
  successful
* The `prior_trunc` option of `identification` was not working
* The `rand_multivariate_student` value of the `proposal_distribution` option
  was not working when used with the
  `tailored_random_block_metropolis_hastings` posterior sampling method
* Perfect foresight simulations of backward models would crash if convergence
  failed with complex-valued residuals
* The diffuse Kalman smoother would crash if `Finf` became singular

Announcement for Dynare 4.6.1 (on 2020-03-13)
=============================================

We are pleased to announce the release of Dynare 4.6.1.

This maintenance release fixes various bugs.

The Windows, macOS and source packages are already available for download at
[the Dynare website](https://www.dynare.org/download/).

All users are strongly encouraged to upgrade.

This release is compatible with MATLAB versions ranging from 7.9 (R2009b) to
9.7 (R2019b), and with GNU Octave versions 5.2.0 (under Windows) and 4.4.1
(under macOS).

Here is a list of the problems identified in version 4.6.0 and that have been
fixed in version 4.6.1:

* Installation on macOS would fail if the GCC compiler was supposed to be
  installed and `www.google.com` was not reachable or blocked
* Dynare++ was missing the `dynare_simul.m` file
* The parameter vector `M_.params` would not be correctly updated after calls
  to `stoch_simul` and `discretionary_policy` if parameters had been modified
  in a steady state file
* The `stoch_simul` command with both the `nograph` and `TeX` options would
  crash
* The `stoch_simul` command with the `noprint` option would crash
* The `prior moments` command would crash if the used parameter vector
  triggered an error code
* In case of problem, the `discretionary_policy` command would crash instead of
  aborting with a proper error code
* Computing of prior/posterior statistics would not work in parallel
* Closing of parallel estimation on GNU/Linux could crash
* The `histval` command would not work in combination with the
  `predetermined_variables` command
* Ramsey optimal policy with multiple instruments would crash if a steady state
  file returned complex values, instead of providing an error message
* The `model_diagnostics` command would not correctly update the parameter
  vector if the latter was set in a steady state file
* The `model_diagnostics` command would ignore the `nocheck` steady state flag


Announcement for Dynare 4.6.0 (on 2020-02-20)
=============================================

We are pleased to announce the release of Dynare 4.6.0.

This major release adds new features and fixes various bugs.

The Windows, macOS and source packages are already available for download at
[the Dynare website](https://www.dynare.org/download/).

All users are strongly encouraged to upgrade.

This release is compatible with MATLAB versions ranging from 7.9 (R2009b) to
9.7 (R2019b), and with GNU Octave versions 5.2.0 (under Windows) and 4.4.1
(under macOS).

Major user-visible changes
--------------------------

 - Stochastic simulations

    - The perturbation method is now available at an arbitrary approximation
      order. In other words, the `order` option of `stoch_simul` accepts an
      arbitrary positive integer (of course, up to some model-specific
      computational limit).

    - New option `filtered_theoretical_moments_grid` of `stoch_simul`, that
      supersedes `hp_ngrid`.

 - Estimation

    - Nonlinear estimation is now also available at an arbitrary approximation
      order. In other words, the `order` option of `estimation` accepts an
      arbitrary positive integer (of course, up to some model-specific
      computational limit).

    - Various improvements to particle filters.

    - It is now possible to estimate models under optimal policy (see below).

    - Variance decomposition of observables now accounts for measurement error.

    - New option `mh_tune_jscale` of `estimation` command for tuning the scale
      parameter of the proposal distribution of the Random Walk Metropolis
      Hastings.

    - Added debugging info when parameters take a `NaN` or `Inf` value.

    - Option `mode_compute=1` is now available under Octave.

 - Perfect foresight and extended path

    - A significant speed improvement should be noted on large models (when
      neither `bytecode` nor `block` option is used). The stacked problem is
      now constructed using a dedicated machine-compiled library that greatly
      speeds up the process (in particular, the time spent in that step can
      become negligible when the `use_dll` option is used).

    - New options `print` and `noprint` of `perfect_foresight_solver` command.

    - Option `stack_solve_algo=2` is now available under Octave.

 - Steady state

    - Option `solve_algo=7` is now available under Octave.

 - Optimal policy

    - The `ramsey_policy` command is now deprecated. It is superseded by
      successive calls to `ramsey_model`, `stoch_simul`, and
      `evaluate_planner_objective` (in this order).

    - It is now possible to estimate a model under optimal policy (either
      Ramsey or discretionary) by running the `estimation` command after either
      `ramsey_model` or `discretionary_policy`. It is however not yet possible
      to estimate parameters that appear in the discount factor of the social
      planner.

    - Discretionary policy returns a more informative error message when the
      objective has nonzero derivatives with respect to some variables.

 - Identification

    - Added minimal system identification check of *Komunjer and Ng (2011)*.

    - Added spectrum identification check of *Qu and Tkachenko (2012)*.

    - Identification is now also available for approximation orders 2 and 3
      with either analytical or numerical parameter derivatives. The relevant
      moments and spectrum are computed from the pruned state space system
      as in *Mutschler (2015)*.

    - All tests (moments, spectrum, minimal system, strength) can be turned
      off.

    - More numerical options can be changed by the user.

    - Improved printing and storage (same folder) of results.

 - Sensitivity analysis

    - New `diffuse_filter` option to the `dynare_sensitivity` command.

    - Arbitrary expressions can now be passed for the interval boundaries in
      `irf_calibration` and `moment_calibration`. ⚠ This breaks the
      previous syntax, requiring that the lower/upper bounds be separated by
      commas.

 - Forecasting and smoothing

    - In `conditional_forecast_paths`, it is no longer required that all
      constrained paths be of the same length. There may now be a different
      number of controlled variables at each period. In that case, the order of
      declaration of endogenous controlled variables and of `controlled_varexo`
      matters: if the second endogenous variable is controlled for less periods
      than the first one, the second `controlled_varexo` isn't set for the last
      periods.

    - New option `parameter_set` to the `calib_smoother` command.

    - ⚠ The results of `conditional_forecast` command is now saved in
      `oo_` (used to be in a file)

 - Shock decomposition

   - Added `fast_realtime` option to real time shock decomposition (deactivated
     by default, runs the smoother only twice: once for the last in-sample and
     once for the last out-of-sample data point).

   - New `diff`, `flip`, `max_nrows`, `plot_init_date` and `plot_end_date`
     options to `plot_shock_decomposition`.

   - New `initial_decomposition_decomposition` command, for computing and
     plotting the decomposition of the effect of smoothed initial conditions of
     state variables.

   - New `squeeze_shock_decomposition` command, for removing decompositions of
     variables that are not of interest.

   - New `with_epilogue` option (common to `shock_decomposition`,
     `realtime_shock_decomposition` and `initial_condition_decomposition`).

   - New `init2shocks` block to attribute initial conditions to shocks.

 - Macro processor

   - New object types: real (supersedes integers), boolean (distinct from
     integers), tuple, user-defined function.

   - New operators: various mathematical functions, set operations on arrays
     (union, intersection, difference, cartesian power and product), type
     checking and conversion.

   - Added support for comprehensions (*e.g.* the set containing the squares of
     all even numbers between 1 and 5 can be constructed with `[ i^2 for i in
     1:5 when mod(i,2) == 0]`).

   - User-defined functions can be declared using the `@#define` operator (*e.g.*
     `@#define f(x) = 2*x^2+3*x+5`).

   - `@#elseif`-clauses are now supported in conditional statements.

   - `@#for` loops can iterate over several variables at the same time (*e.g.*
     `@#for (i,j) in X`, where `X` is an array containing tuples of size 2).

   - Added the possibility to exclude some elements when iterating over `@#for`
     loops (*e.g.* `@#for i in 1:5 when mod(i,2) == 0` iterates over all even
     numbers between 1 and 5).

   - A `defined()` function allows testing whether macro variables have been
     defined.

   - Empty arrays (with the `[]` syntax) are now possible.

   - Arrays of arrays are now supported.

   - New macro directives `@#echomacrovars` and `@#echomacrovars(save)` for
     displaying or saving the values of all macro-variables.

   - Inline comments are now supported.

   - ⚠ All division operations are now done with doubles (as opposed to
     integers). To achieve the old functionality, use the new `floor` operator.

   - ⚠ Colon syntax used to require braces around it to create an array
     (*e.g.* `[1:3]` would create `[1,2,3]`). Now this is not necessary (`1:3`
     creates `[1,2,3]` while `[1:3]` would create `[[1,2,3]]`).

   - ⚠ Previously, printing a boolean would print `1` or `0`. Now, it
     prints `true` or `false`. To achieve the old functionality, you must cast
     it to a real, *e.g.* `@{(real)(1!=0)}`.

 - LaTeX output

    - New command `write_latex_steady_state_model`.

    - New option `planner_discount_latex_name` of `ramsey_model` and
      `discretionary_policy`.

    - New command `model_local_variable` command for assigning a LaTeX name to
      model-local variables.

    - The `write_latex_static_model` and `write_latex_original_model` commands
      now support the `write_equation_tags` option.

 - Compilation of the model (`use_dll` option) made easier and faster

    - Under Windows, it is no longer necessary to manually install the
      compiler, since the latter is now shipped by the Dynare installer.

    - Under macOS, the Dynare installer now automatically downloads and
      installs the compiler.

    - It is no longer necessary to configure MATLAB to let it know where the
      compiler is, since the compilation is now done by the preprocessor.

    - The compilation phase is now faster on large models (this has been
      achieved by disabling a few time-consuming and not-so-useful optimization
      passes otherwise done by the compiler).

    - New `compilation_setup` block for specifying a custom compiler or custom
      compilation flags.

 - Model, variables and parameters declaration

    - New syntax to declare model variables and parameters on-the-fly in the
      `model` block. To do this, simply follow the symbol name with a vertical
      line (`|`, pipe character) and either an `e`, an `x`, or a `p`. For
      example, to declare a parameter named `alpha` in the model block, you
      could write `alpha|p` directly in an equation where it appears.
      Similarly, to declare an endogenous variable `c` in the model block you
      could write `c|e`.

    - New syntax to declare model variable and parameters on-the-fly in
      equation tags. In the tag, simply state the type of variable to be
      declared (`endogenous`, `exogenous`, or `parameter` followed by an equal
      sign and the variable name in single quotes. Hence, to declare a variable
      `c` as endogenous in an equation tag, you can type `[endogenous='c']`.

    - New `epilogue` block for computing output variables of interest that may
      not be necessarily defined in the model (*e.g.* various kinds of
      real/nominal shares or relative prices, or annualized variables out of a
      quarterly model).

 - Command-line options

    - Added the possibility to declare Dynare command-line options in the `.mod`
      file.

    - New option `nopreprocessoroutput` to disable printing of messages from
      the preprocessor.

    - It is now possible to assign an arbitrary macro-expression to a
      macro-variable defined on the command-line, using the `-D` syntax.

    - New  option `linemacro` to revert to the old format of the
      macro-processed file (see below).

 - Preprocessor outputs and inputs

    - Added JSON output to the preprocessor. A representation of the model file
      and the whole content of the `.mod` file is saved in `.json` files.
      These JSON files can be easily parsed from any language (C++, Fortran,
      Python, Julia, MATLAB, Octave…). This new feature opens the possibility to
      develop alternative back-ends for the Dynare language.

    - ⚠ Most files generated by the preprocessor are now grouped under
      two subdirectories. Assuming your file is `FILENAME.mod`, then M-files
      and MEX-files will be under `+FILENAME/`, while other output (JSON,
      LaTeX, source code for the MEX files) will be under `FILENAME/`.

    - The macro-generated output is now more readable (no more line numbers and
      empty lines). The old behaviour can be restored using the `linemacro`
      option (see above).

    - Ability to call the preprocessor by passing the `.mod` file as a string
      argument from the macOS or GNU/Linux command line.

 - dseries classes

    - New functionalities and efficiency improvements.

    - Complete rewrite using the new `classdef` syntax and exploiting in place
      modifications when possible.

    - Integration of the `dates` classes within `dseries`.

 - Reporting classes

    - Automatically create titlepage with page numbers/page titles.

    - Allow for the removal of headers and footers from a given page.

    - Allow user to set page number.

    - Split up report output. Create new files for the preamble, the body of
      the report, and each individual page of the report.

    - The classes have been converted to the new `classdef` syntax.

 - Misc

    - External functions can be located in MATLAB/Octave namespaces.

    - Improvements to the balanced growth path test that is performed after
      Dynare has detrended the model (given the trends on variables declared by
      the user): the default tolerance has been raised, and a different value
      can be set with new option `balanced_growth_test_tol` to the `model`
      block; as a consequence, failing the test is now an error again.

    - New collection of MATLAB/Octave utilities to retrieve and alter objects:
      `get_irf`, `get_mean`, `get_shock_stderr_by_name`, `get_smooth`,
      `get_update`, `set_shock_stderr_value`.

    - ⚠ Previously, when some MEX files were missing, Dynare would
      automatically fall back to slower M-file functional alternative; this is
      no longer the case. It is however still possible to manually add these
      alternatives in the MATLAB/Octave path (they are located under
      `matlab/missing/mex`; this only applies to the `mjdgges`, `gensylv`,
      `A_times_B_kronecker_C`, `sparse_hessian_times_B_kronecker_C` and
      `local_state_space_iteration_2` DLLs).


Since there are a few backward-incompatible changes in this release, users may
want to have a look at the [upgrade
guide](https://git.dynare.org/Dynare/dynare/-/wikis/BreakingFeaturesIn4.6) to
adapt their existing codes.


Bugs that were present in 4.5.7 and that are fixed in 4.6.0
-----------------------------------------------------------

* Estimation: the check for stochastic singularity erroneously would only take
  estimated measurement error into account.
* Estimation: if the Hessian at the mode was not positive definite, the Laplace
  approximation returned a complex number, but only displayed the real-valued
  part.
* Conditional Forecasting: using one period only would result in a crash.
* First-order approximation was not working with purely forward-looking models.
* The preprocessor would not allow for inline comments including macro
  statements.
* Using the `STEADY_STATE()` operator on exogenous variables would lead to
  crashes in stochastic simulations.
* `moment_calibration`: for autocorrelation functions, the x-axis labeling had
  the wrong order.
* `plot_identification`: placement of white dots indicating infinite values was
  incorrect
* Automatic detrending would sometime refuse to detrend model despite the user
  having given correct trends.
* Using `use_dll` + `fast` options would not always recompile the model when
  the equations were changed.
* Under certain circumstances, the combination of `bytecode` and
  `stack_solve_algo=1` options could lead to crashes or wrong results.


References
----------

 - Komunjer, I. and S. Ng (2011), “[Dynamic Identification of Dynamic
   Stochastic General Equilibrium
   Models](https://www.onlinelibrary.wiley.com/doi/abs/10.3982/ECTA8916),”
   *Econometrica*, 79(6), 1995–2032

 - Qu, Z. and D. Tkachenko (2012), “[Identification and frequency domain
   quasi‐maximum likelihood estimation of linearized dynamic stochastic
   general equilibrium
   models](https://onlinelibrary.wiley.com/doi/abs/10.3982/QE126),”
   *Quantitative Economics*, 3(1), 95–132

 - Mutschler, W. (2015), “[Identification of DSGE models—The effect of
   higher-order approximation and
   pruning](https://www.sciencedirect.com/science/article/pii/S0165188915000731),”
   *Journal of Economic Dynamics and Control*, 56, 34–54


Announcement for Dynare 4.5.7 (on 2019-02-06)
=============================================

We are pleased to announce the release of Dynare 4.5.7.

This is a bugfix release.

The Windows packages are already available for download at: <http://www.dynare.org/download/dynare-stable>.

The Mac and GNU/Linux packages (for Debian and Ubuntu LTS) should follow soon.

This release is compatible with MATLAB versions 7.5 (R2007b) to 9.4 (R2018b)
and with GNU Octave versions 4.4.1.

Here is a list of the problems identified in version 4.5.6 and that have been
fixed in version 4.5.7:

 - The mex-file conducting the QZ decomposition erroneously applied
   the `qz_criterium` to the square absolute value of eigenvalues
   instead of the absolute value itself (as done in mjdgges.m and the
   AIM solver).

 - In pathological cases, `mode_compute=5` (`newrat`) might enter an
   infinite loop.

 - `discretionary_policy` might erroneously state that the derivatives
   of the objective function are non-zero if there are NaN present.

 - Dynare++, when conducting the QZ decomposition, erroneously applied
   the `qz_criterium` to the square absolute value of eigenvalues
   instead of the absolute value itself.

 - Dynare++: IRFs were incorrectly computed.

 - `dynare_sensitivity` did not display the figures of
   `irf_calibration`, it only stored them on the disk.

 - Scatter plots generated by `dynare_sensitivity` did not correctly
   display LaTeX names.

 - Parameter updating via steady state files did not correctly work in
   case of using `[static]`/`[dynamic]` equation tags.

 - Memory leaks in `k_order_pert` (used by higher order stochastic
   simulations) could lead to crashes.

 - Predetermined variables were not properly set when used in model
   local variables.

 - Posterior moment computation did not correctly update the
   covariance matrix of exogenous shocks during posterior sampling.

 - Dynare was crashing with a cryptic message if a non estimated
   parameter was initialized in the `estimated_params_init` block.

 - The `forecast` command crashed if the model was declared as linear
   and contained deterministic exogenous variables.

 - Block decomposition is broken when used in conjunction with
   `varexo_det`.

 - The model was not correctly specified when `identification` was run
   without another stochastic command in the `.mod` file
   (*e.g.* `estimation`, `stoch_simul`, etc.).

 - Realtime annualized shock decompositions added the wrong steady state
   value.

 - `mh_recover` option crashed when using slice sampler.

 - x-axis values in plots of moment restrictions were wrong for
   autocovariances.



Announcement for Dynare 4.5.6 (on 2018-07-25)
=============================================

We are pleased to announce the release of Dynare 4.5.6.

This is a bugfix release.

The Windows packages are already available for download at:
<http://www.dynare.org/download/dynare-stable>.

The Mac and GNU/Linux packages (for Debian and Ubuntu LTS) should follow soon.

This release is compatible with MATLAB versions 7.5 (R2007b) to 9.4 (R2018a)
and with GNU Octave versions 4.4.

Here is a list of the problems identified in version 4.5.5 and that have been
fixed in version 4.5.6:

 - TaRB sampler: incorrect last posterior was returned if the last draw was
   rejected.

 - Fixed online particle filter by drawing initial conditions in the prior
   distribution.

 - Fixed evaluation of the likelihood in non linear / particle filters.

 - Added missing documented `montecarlo` option in Gaussian Filter and
   Nonlinear Kalman Filter.

 - Added back a flag to deal with errors on Cholesky decomposition in the
   Conditional Particle Filter.

 - Macroprocessor `length()` operator was returning 1 when applied to a
   string. Macroprocessor now raises an error when `length()` operator is
   called on an integer and return the number of characters when applied to a
   string.

 - `mode_compute=8`: the error code during mode-finding was not correctly
   handled, resulting in crashes.

 - Identification was not correctly displaying a message for collinear parameters
   if there was no unidentified parameter present.



Announcement for Dynare 4.5.5 (on 2018-06-08)
=============================================

We are pleased to announce the release of Dynare 4.5.5.

This is a bugfix release.

The Windows packages are already available for download at:
<http://www.dynare.org/download/dynare-stable>.

The Mac and GNU/Linux packages (for Debian and Ubuntu LTS) should follow soon.

This release is compatible with MATLAB versions 7.5 (R2007b) to 9.4 (R2018a)
and with GNU Octave versions 4.2.

Here is a list of the problems identified in version 4.5.4 and that have been
fixed in version 4.5.5:

 - Identification was crashing during prior sampling if `ar` was initially too
   low.

 - The `align` method on `dseries` did not return a functional second `dseries`
   output.

 - Predetermined variables were not properly set when used in model local
   variables.

 - `perfect_foresight_solver` with option `stack_solve_algo=7` was not working
   correctly when an exogenous variable has a lag greater than 1.

 - `identification` with `prior_mc` option would crash if the number of moments
   with non-zero derivative is smaller than the number of parameters.

 - Calling several times `normcdf` or `normpdf` with the same arguments in a
   model with block decomposition (but not bytecode) was leading to incorrect
   results.



Announcement for Dynare 4.5.4 (on 2018-01-29)
=============================================

We are pleased to announce the release of Dynare 4.5.4.

This is a bugfix release.

The Windows packages are already available for download at:
<http://www.dynare.org/download/dynare-stable>.

The Mac and GNU/Linux packages (for Debian and Ubuntu LTS) should follow soon.

This release is compatible with MATLAB versions 7.5 (R2007b) to 9.3 (R2017b)
and with GNU Octave versions 4.2.

Here is a list of the problems identified in version 4.5.3 and that have been
fixed in version 4.5.4:

 - The `type` option of `plot_shock_decomposition` was always set to `qoq` regardless of what is specified.

 - Bug in GSA when no parameter was detected below pvalue threshold.

 - Various bug fixes in shock decompositions.

 - Bug in reading in macro arrays passed on `dynare` command line via the `-D` option.

 - Estimation with missing values was crashing if the `prefilter` option was used.

 - Added a workaround for a difference in behaviour between Octave and MATLAB regarding the creation
   of function handles for functions that do not exist in the path. With Octave 4.2.1, steady state
   files did not work if no auxiliary variables were created.

 - The `stoch_simul` command was crashing with a cryptic message if option `order=3` was used without
   setting `k_order_solver`.

 - In cases where the prior bounds are infinite and the mode is estimated at exactly 0, no `mode_check`
   graphs were displayed.

 - Parallel execution of MCMC was broken in models without auxiliary variables.

 - Reading data with column names from Excel might crash.

 - The multivariate Kalman smoother was crashing in case of missing data in the observations and
   `Finf` became singular.

 - The `plot_shock_decomposition` command ignored various user-defined options like `fig_name`,
   `use_shock_groups` or `interactive` and instead used the default options.

 - Nested `@#ifdef` and `@#ifndef` statements don’t work in the macroprocessor.



Announcement for Dynare 4.5.3 (on 2017-10-19)
=============================================

We are pleased to announce the release of Dynare 4.5.3.

This is a bugfix release. It comes less than 24 hours after the previous release,
because version 4.5.2 was affected by a critical bug for MATLAB older than R2016b.

The Windows packages are already available for download at:
<http://www.dynare.org/download/dynare-stable>.

The Mac and GNU/Linux packages (for Debian and Ubuntu LTS) should follow soon.

This release is compatible with MATLAB versions 7.5 (R2007b) to 9.3 (R2017b)
and with GNU Octave versions 4.2.

Here is a list of the problems identified in version 4.5.2 and that have been
fixed in version 4.5.3:


 - `isfile` routine was failing with MATLAB older than R2016b. This bug did not
   affect Octave.



Announcement for Dynare 4.5.2 (on 2017-10-19)
=============================================

We are pleased to announce the release of Dynare 4.5.2.

This is a bugfix release.

The Windows packages are already available for download at:
<http://www.dynare.org/download/dynare-stable>.

The Mac and GNU/Linux packages (for Debian and Ubuntu LTS) should follow soon.

This release is compatible with MATLAB versions 7.5 (R2007b) to 9.3 (R2017b)
and with GNU Octave versions 4.2.

Here is a list of the problems identified in version 4.5.1 and that have been
fixed in version 4.5.2:


 - Fixed bug in perfect foresight solver:

   + If expected shocks were declared after the terminal period, as specified
   by the `periods` option, Dynare was crashing.

   + Models declared with the `linear` option were crashing if exogenous
   variables were present with a lead or lag.

 - After ML or Bayesian estimation when the smoother option or `mh_replic=0`
   were not specified, not all smoothed measurement errors were displayed.

 - Fixed error in reference manual about the `conditional_forecasts` command.

 - Fixed smoother behaviour, provide informative error instead of crashing when
   model cannot be solved.

 - The `nopathchange` preprocessor option was always triggered, regardless of
   whether it was passed or not.

 - When `ramsey_policy` is used, allow state variables to be set in `histval`
   block.

 - `histval` erroneously accepted leads, leading to cryptic crashes.

 - The prior MC draws from previous runs were not deleted, potentially
   resulting in loading stale files.

 - `estim_params_` was being declared `global` more than once.

 - Fixed crashes happening when simulating linear models with order>1.

 - Make empirical moments independent of `simul_replic`, as stated in the
   reference manual, by outputting moments computed with the first simulated
   sample.

 - The `prior_function` required a preceding `estimation`-command to properly
   set up the prior.

 - If the mode for a parameter was at exactly 0, `mode_check` was crashing.

 - Fixed `get_posterior_parameters`-routine which should not do more than
   getting parameters. As a consequense, the `shock_decomposition`-command
   did not correctly set the `parameter_set` for use in subsequent function
   calls if shocks are correlated or measurement error is present.

 - Fixed bug in Ramsey problem with constraints both on a policy instrument and
   another variable. Note that the constraint on a variable that is not an
   instrument of the Ramsey problem must be written with an equation tag in the
   model block.

 - Fixed bug in Ramsey problem with constraints on policy instrument.

 - Fixed crash with optimizer 5 when not used with DSGE model at order 1.

 - Fixed mex file used for third order approximation (was crashing on
   MATLAB/Windows 7).



Announcement for Dynare 4.5.1 (on 2017-08-24)
=============================================

We are pleased to announce the release of Dynare 4.5.1.

This is a bugfix release.

The Windows packages are already available for download at:
<http://www.dynare.org/download/dynare-stable>.

The Mac and GNU/Linux packages (for Debian and Ubuntu LTS) should follow soon.

This release is compatible with MATLAB versions 7.5 (R2007b) to 9.2 (R2017a)
and with GNU Octave versions 4.2.

Here is a list of the problems identified in version 4.5.0 and that have been
fixed in version 4.5.1:


 - Fixed out of memory issue with simpsa optimization algorithm.

 - Added missing plots for measurement errors with `generate_trace_plot`
   command.

 - Posterior moments after MCMC for very big models were not correctly computed
   and their plotting might crash Dynare.

 - Results of the posterior conditional variance decomposition after MCMC were
   not correctly computed.

 - Options `use_shock_groups` and `colormap` of the `shock_decomposition`
   command were not working.

 - Added a clean error message if sensitivity toolbox is used with recursive
   estimation.

 - Computation of posterior filtered variables was crashing in models with only
   one variable.

 - Fixed various typos and errors in the reference manual.



Announcement for Dynare 4.5.0 (on 2017-06-11)
=============================================

We are pleased to announce the release of Dynare 4.5.0.

This major release adds new features and fixes various bugs.

The Windows packages are already available for download at:
<http://www.dynare.org/download/dynare-stable>.

The Mac and Debian/Ubuntu packages should follow soon.

All users are strongly encouraged to upgrade.

This release is compatible with MATLAB versions ranging from 7.5 (R2007b) to
9.2 (R2017a) and with GNU Octave version 4.2.

Here is the list of major user-visible changes:


 - Ramsey policy

   + Added command `ramsey_model` that builds the expanded model with
     FOC conditions for the planner’s problem but doesn’t perform any
     computation. Usefull to compute Ramsey policy in a perfect
     foresight model,

   + `ramsey_policy` accepts multipliers in its variable list and
     displays results for them.


 - Perfect foresight models

   + New commands `perfect_foresight_setup` (for preparing the
     simulation) and `perfect_foresight_solver` (for computing it). The
     old `simul` command still exist and is now an alias for
     `perfect_foresight_setup` + `perfect_foresight_solver`. It is no
     longer possible to manipulate by hand the contents of
     `oo_.exo_simul` when using `simul`. People who want to do
     it must first call `perfect_foresight_setup`, then do the
     manipulations, then call `perfect_foresight_solver`,

   + By default, the perfect foresight solver will try a homotopy
     method if it fails to converge at the first try. The old behavior
     can be restored with the `no_homotopy` option,

   + New option `stack_solve_algo=7` that allows specifying a
     `solve_algo` solver for solving the model,

   + New option `solve_algo` that allows specifying a solver for
     solving the model when using `stack_solve_algo=7`,

   + New option `lmmcp` that solves the model via a Levenberg-Marquardt
     mixed complementarity problem (LMMCP) solver,

   + New option `robust_lin_solve` that triggers the use of a robust
     linear solver for the default `solve_algo=4`,

   + New options `tolf` and `tolx` to control termination criteria of
     solvers,

   + New option `endogenous_terminal_period` to `simul`,

   + Added the possibility to set the initial condition of the
     (stochastic) extended path simulations with the histval block.


 - Optimal simple rules

   + Saves the optimal value of parameters to `oo_.osr.optim_params`,

   + New block `osr_params_bounds` allows specifying bounds for the
     estimated parameters,

   + New option `opt_algo` allows selecting different optimizers while
     the new option `optim` allows specifying the optimizer options,

   + The `osr` command now saves the names, bounds, and indices for the
     estimated parameters as well as the indices and weights of the
     variables entering the objective function into `M_.osr`.


 - Forecasts and Smoothing

   + The smoother and forecasts take uncertainty about trends and means
     into account,

   + Forecasts accounting for measurement error are now saved in fields
     of the form `HPDinf_ME` and `HPDsup_ME`,

   + New fields `oo_.Smoother.Trend` and `oo_.Smoother.Constant` that
     save the trend and constant parts of the smoothed variables,

   + new field `oo_.Smoother.TrendCoeffs` that stores the trend
     coefficients.

   + Rolling window forecasts allowed in `estimation` command by
     passing a vector to `first_obs`,

   + The `calib_smoother` command now accepts the `loglinear`,
     `prefilter`, `first_obs` and `filter_decomposition` options.


 - Estimation

   + New options: `logdata`, `consider_all_endogenous`,
     `consider_only_observed`, `posterior_max_subsample_draws`,
     `mh_conf_sig`, `diffuse_kalman_tol`, `dirname`, `nodecomposition`

   + `load_mh_file` and `mh_recover` now try to load chain’s proposal density,

   + New option `load_results_after_load_mh` that allows loading some
     posterior results from a previous run if no new MCMC draws are
     added,

   + New option `posterior_nograph` that suppresses the generation of
     graphs associated with Bayesian IRFs, posterior smoothed objects,
     and posterior forecasts,

   + Saves the posterior density at the mode in
     `oo_.posterior.optimization.log_density`,

   + The `filter_covariance` option now also works with posterior
     sampling like Metropolis-Hastings,

   + New option `no_posterior_kernel_density` to suppress computation
     of kernel density of posterior objects,

   + Recursive estimation and forecasting now provides the individual
     `oo_` structures for each sample in `oo_recursive_`,

   + The `trace_plot` command can now plot the posterior density,

   + New command `generate_trace_plots` allows generating all trace
     plots for one chain,

   + New commands `prior_function` and `posterior_function` that
     execute a user-defined function on parameter draws from the
     prior/posterior distribution,

   + New option `huge_number` for replacement of infinite bounds with
     large number during `mode_compute`,

   + New option `posterior_sampling_method` allows selecting the new
     posterior sampling options:
     `tailored_random_block_metropolis_hastings` (Tailored randomized
     block (TaRB) Metropolis-Hastings), `slice` (Slice sampler),
     `independent_metropolis_hastings` (Independent
     Metropolis-Hastings),

   + New option `posterior_sampler_options` that allow controlling the
     options of the `posterior_sampling_method`, its `scale_file`-option
     pair allows loading the `_mh_scale.mat`-file storing the tuned
     scale factor from a previous run of `mode_compute=6`,

   + New option `raftery_lewis_diagnostics` that computes *Raftery and Lewis
     (1992)* convergence diagnostics,

   + New option `fast_kalman_filter` that provides fast Kalman filter
     using Chandrasekhar recursions as described in *Ed Herbst (2015)*,

   + The `dsge_var` option now saves results at the posterior mode into
     `oo_.dsge_var`,

   + New option `smoothed_state_uncertainty` to provide the uncertainty
     estimate for the smoothed state estimate from the Kalman smoother,

   + New prior density: generalized Weibull distribution,

   + Option `mh_recover` now allows continuing a crashed chain at the
     last save mh-file,

   + New option `nonlinear_filter_initialization` for the
     `estimation` command. Controls the initial covariance matrix
     of the state variables in nonlinear filters.

   + The `conditional_variance_decomposition` option now displays
     output and stores it as a LaTeX-table when the `TeX` option is
     invoked,

   + The `use_calibration` to `estimated_params_init` now also works
     with ML,

   + Improved initial estimation checks.


 - Steady state

   + The default solver for finding the steady state is now a
     trust-region solver (can be triggered explicitly with option
     `solve_algo=4`),

   + New options `tolf` and `tolx` to control termination criteria of
     solver,

   + The debugging mode now provides the termination values in steady
     state finding.


 - Stochastic simulations

   + New options `nodecomposition`,

   + New option `bandpass_filter` to compute bandpass-filtered
     theoretical and simulated moments,

   + New option `one_sided_hp_filter` to compute one-sided HP-filtered
     simulated moments,

   + `stoch_simul` displays a simulated variance decomposition when
     simulated moments are requested,

   + `stoch_simul` saves skewness and kurtosis into respective fields
     of `oo_` when simulated moments have been requested,

   + `stoch_simul` saves the unconditional variance decomposition in
     `oo_.variance_decomposition`,

   + New option `dr_display_tol` that governs omission of small terms
     in display of decision rules,

   + The `stoch_simul` command now prints the displayed tables as LaTeX
     code when the new `TeX` option is enabled,

   + The `loglinear` option now works with lagged and leaded exogenous
     variables like news shocks,

   + New option `spectral_density` that allows displaying the spectral
     density of (filtered) endogenous variables,

   + New option `contemporaneous_correlation` that allows saving
     contemporaneous correlations in addition to the covariances.


 - Identification

   + New options `diffuse_filter` and `prior_trunc`,

   + The `identification` command now supports correlations via
     simulated moments,


 - Sensitivity analysis

   + New blocks `irf_calibration` and `moment_calibration`,

   + Outputs LaTeX tables if the new `TeX` option is used,

   + New option `relative_irf` to `irf_calibration` block.


 - Conditional forecast

   + Command `conditional_forecast` now takes into account `histval`
     block if present.


 - Shock decomposition

   + New option `colormap` to `shocks_decomposition` for controlling
     the color map used in the shocks decomposition graphs,

   + `shocks_decomposition` now accepts the `nograph` option,

   + New command `realtime_shock_decomposition` that for each period `T= [presample,...,nobs]`
     allows computing the:

     * realtime historical shock decomposition `Y(t|T)`, *i.e.* without observing data in `[T+1,...,nobs]`

     * forecast shock decomposition `Y(T+k|T)`

     * realtime conditional shock decomposition `Y(T+k|T+k)-Y(T+k|T)`

   + New block `shock_groups` that allows grouping shocks for the
     `shock_decomposition` and `realtime_shock_decomposition` commands,

   + New command `plot_shock_decomposition` that allows plotting the
     results from `shock_decomposition` and
     `realtime_shock_decomposition` for different vintages and shock
     groupings.


 - Macroprocessor

   + Can now pass a macro-variable to the `@#include` macro directive,

   + New preprocessor flag `-I`, macro directive `@#includepath`, and
     dynare config file block `[paths]` to pass a search path to the
     macroprocessor to be used for file inclusion via `@#include`.


 - Command line

   + New option `onlyclearglobals` (do not clear JIT compiled functions
     with recent versions of MATLAB),

   + New option `minimal_workspace` to use fewer variables in the
     current workspace,

   + New option `params_derivs_order` allows limiting the order of the
     derivatives with respect to the parameters that are calculated by
     the preprocessor,

   + New command line option `mingw` to support the MinGW-w64 C/C++
     Compiler from TDM-GCC for `use_dll`.


 - dates/dseries/reporting classes

   + New methods `abs`, `cumprod` and `chain`,

   + New option `tableRowIndent` to `addTable`,

   + Reporting system revamped and made more efficient, dependency on
     matlab2tikz has been dropped.


 - Optimization algorithms

   + `mode_compute=2` Uses the simulated annealing as described by
     *Corana et al. (1987)*,

   + `mode_compute=101` Uses SOLVEOPT as described by *Kuntsevich and
     Kappel (1997)*,

   + `mode_compute=102` Uses `simulannealbnd` from MATLAB’s Global
     Optimization Toolbox (if available),

   + New option `silent_optimizer` to shut off output from mode
     computing/optimization,

   + New options `verbosity` and `SaveFiles` to control output and
     saving of files during mode computing/optimization.


 - LaTeX output

   + New command `write_latex_original_model`,

   + New option `write_equation_tags` to `write_latex_dynamic_model`
     that allows printing the specified equation tags to the generate
     LaTeX code,

   + New command `write_latex_parameter_table` that writes the names and
     values of model parameters to a LaTeX table,

   + New command `write_latex_prior_table` that writes the descriptive
     statistics about the prior distribution to a LaTeX table,

   + New command `collect_latex_files` that creates one compilable LaTeX
     file containing all TeX-output.


 - Misc.

   + Provides 64bit preprocessor,

   + Introduces new path management to avoid conflicts with other
     toolboxes,

   + Full compatibility with MATLAB 2014b’s new graphic interface,

   + When using `model(linear)`, Dynare automatically checks
     whether the model is truly linear,

   + `usedll`, the `msvc` option now supports `normcdf`, `acosh`,
     `asinh`, and `atanh`,

   + New parallel option `NumberOfThreadsPerJob` for Windows nodes that
     sets the number of threads assigned to each remote MATLAB/Octave
     run,

   + Improved numerical performance of
     `schur_statespace_transformation` for very large models,

   + The `all_values_required` option now also works with `histval`,

   + Add missing `horizon` option to `ms_forecast`,

   + BVAR now saves the marginal data density in
     `oo_.bvar.log_marginal_data_density` and stores prior and
     posterior information in `oo_.bvar.prior` and
     `oo_.bvar.posterior`.



Bugs and problems identified in version 4.4.3 and that have been fixed in version 4.5.0:


 - BVAR models

   + `bvar_irf` could display IRFs in an unreadable way when they moved from
     negative to positive values,

   + In contrast to what is stated in the documentation, the confidence interval
     size `conf_sig` was 0.6 by default instead of 0.9.


 - Conditional forecasts

   + The `conditional_forecast` command produced wrong results in calibrated
     models when used at initial values outside of the steady state (given with
     `initval`),

   + The `plot_conditional_forecast` option could produce unreadable figures if
     the areas overlap,

   + The `conditional_forecast` command after MLE crashed,

   + In contrast to what is stated in the manual, the confidence interval size
     `conf_sig` was 0.6 by default instead of 0.8.

   + Conditional forecasts were wrong when the declaration of endogenous
     variables was not preceeding the declaration of the exogenous
     variables and parameters.


 - Discretionary policy

   + Dynare allowed running models where the number of instruments did not match
     the number of omitted equations,

   + Dynare could crash in some cases when trying to display the solution,

   + Parameter dependence embedded via a `steady_state` was not taken into
     account, typically resulting in crashes.

 - dseries class

   + When subtracting a dseries object from a number, the number was instead
     subtracted from the dseries object.


 - DSGE-VAR models

   + Dynare crashed when estimation encountered non-finite values in the Jacobian
     at the steady state,

   + The presence of a constant was not considered for degrees of freedom
     computation of the Gamma function used during the posterior computation; due
     to only affecting the constant term, results should be be unaffected, except
     for model_comparison when comparing models with and without.


 - Estimation command

   + In contrast to what was stated in the manual, the confidence interval size
     `conf_sig` for `forecast` without MCMC was 0.6 by default instead of 0.9,

   + Calling estimation after identification could lead to crashes,

   + When using recursive estimation/forecasting and setting some elements of
     `nobs` to be larger than the number of observations T in the data,
     `oo_recursive_` contained additional cell entries that simply repeated the
     results obtained for `oo_recursive_T`,

   + Computation of Bayesian smoother could crash for larger models when
     requesting `forecast` or `filtered_variables`,

   + Geweke convergence diagnostics were not computed on the full MCMC chain when
     the `load_mh_file` option was used,

   + The Geweke convergence diagnostics always used the default `taper_steps` and
   `geweke_interval`,

   + Bayesian IRFs (`bayesian_irfs` option) could be displayed in an unreadable
     way when they move from negative to positive values,

   + If `bayesian_irfs` was requested when `mh_replic` was too low to compute
     HPDIs, plotting was crashing,

   + The x-axis value in `oo_.prior_density` for the standard deviation and
     correlation of measurement errors was written into a field
     `mearsurement_errors_*` instead of `measurement_errors_*`,

   + Using a user-defined `mode_compute` crashed estimation,

   + Option `mode_compute=10` did not work with infinite prior bounds,

   + The posterior variances and covariances computed by `moments_varendo` were
     wrong for very large models due to a matrix erroneously being filled up with
     zeros,

   + Using the `forecast` option with `loglinear` erroneously added the unlogged
     steady state,

   + When using the `loglinear` option the check for the presence of a constant
     was erroneously based on the unlogged steady state,

   + Estimation of `observation_trends` was broken as the trends specified as a
     function of deep parameters were not correctly updated during estimation,

   + When using `analytic_derivation`, the parameter values were not set before
     testing whether the steady state file changes parameter values, leading to
     subsequent crashes,

   + If the steady state of an initial parameterization did not solve, the
     observation equation could erroneously feature no constant when the
     `use_calibration` option was used,

   + When computing posterior moments, Dynare falsely displayed that moment
     computations are skipped, although the computation was performed correctly,

   + If `conditional_variance_decomposition` was requested, although all
     variables contain unit roots, Dynare crashed instead of providing an error
     message,

   + Computation of the posterior parameter distribution was erroneously based
     on more draws than specified (there was one additional draw for every Markov
     chain),

   + The estimation option `lyapunov=fixed_point` was broken,

   + Computation of `filtered_vars` with only one requested step crashed Dynare,

   + Option `kalman_algo=3` was broken with non-diagonal measurement error,

   + When using the diffuse Kalman filter with missing observations, an additive
     factor log(2π) was missing in the last iteration step,

   + Passing of the `MaxFunEvals` and `InitialSimplexSize` options to
     `mode_compute=8` was broken,

   + Bayesian forecasts contained initial conditions and had the wrong length in
     both plots and stored variables,

   + Filtered variables obtained with `mh_replic=0`, ML, or
     `calibrated_smoother` were padded with zeros at the beginning and end and
     had the wrong length in stored variables,

   + Computation of smoothed measurement errors in Bayesian estimation was broken,

   + The `selected_variables_only` option (`mh_replic=0`, ML, or
     `calibrated_smoother`) returned wrong results for smoothed, updated, and
     filtered variables,

   + Combining the `selected_variables_only` option with forecasts obtained
     using `mh_replic=0`, ML, or `calibrated_smoother` leaded to crashes,

   + `oo_.UpdatedVariables` was only filled when the `filtered_vars` option was specified,

   + When using Bayesian estimation with `filtered_vars`, but without
     `smoother`, then `oo_.FilteredVariables` erroneously also contained filtered
     variables at the posterior mean as with `mh_replic=0`,

   + Running an MCMC a second time in the same folder with a different number of
     iterations could result in crashes due to the loading of stale files,

   + Results displayed after Bayesian estimation when not specifying
     the `smoother` option were based on the parameters at the mode
     from mode finding instead of the mean parameters from the
     posterior draws. This affected the smoother results displayed, but
     also calls to subsequent command relying on the parameters stored
     in `M_.params` like `stoch_simul`,

   + The content of `oo_.posterior_std` after Bayesian estimation was based on
     the standard deviation at the posterior mode, not the one from the MCMC, this
     was not consistent with the reference manual,

   + When the initialization of an MCMC run failed, the metropolis.log file was
     locked, requiring a restart of MATLAB to restart estimation,

   + If the posterior mode was right at the corner of the prior bounds, the
     initialization of the MCMC erroneously crashed,

   + If the number of dropped draws via `mh_drop` coincided with the number of
     draws in a `_mh`-file, `oo_.posterior.metropolis.mean` and
     `oo_.posterior.metropolis.Variance` were NaN.


 - Estimation and calibrated smoother

   + When using `observation_trends` with the `prefilter` option, the mean shift
     due to the trend was not accounted for,

   + When using `first_obs`>1, the higher trend starting point of
     `observation_trends` was not taken into account, leading, among other things,
     to problems in recursive forecasting,

   + The diffuse Kalman smoother was crashing if the forecast error variance
     matrix becomes singular,

   + The multivariate Kalman smoother provided incorrect state estimates when
     all data for one observation are missing,

   + The multivariate diffuse Kalman smoother provided incorrect state estimates
     when the `Finf` matrix becomes singular,

   + The univariate diffuse Kalman filter was crashing if the initial covariance
     matrix of the nonstationary state vector is singular,


 - Forecats

   + In contrast to what is stated in the manual, the confidence interval size
     `conf_sig` was 0.6 by default instead of 0.9.

   + Forecasting with exogenous deterministic variables provided wrong decision
     rules, yielding wrong forecasts.

   + Forecasting with exogenous deterministic variables crashed when the
     `periods` option was not explicitly specified,

   + Option `forecast` when used with `initval` was using the initial values in
     the `initval` block and not the steady state computed from these initial
     values as the starting point of forecasts.


 - Global Sensitivity Analysis

   + Sensitivity with ML estimation could result in crashes,

   + Option `mc` must be forced if `neighborhood_width` is used,

   + Fixed dimension of `stock_logpo` and `stock_ys`,

   + Incomplete variable initialization could lead to crashes with `prior_range=1`.


 - Indentification

   + Identification did not correctly pass the `lik_init` option,
     requiring the manual setting of `options_.diffuse_filter=1` in
     case of unit roots,

   + Testing identification of standard deviations as the only
     parameters to be estimated with ML leaded to crashes,

   + Automatic increase of the lag number for autocovariances when the
     number of parameters is bigger than the number of non-zero moments
     was broken,

   + When using ML, the asymptotic Hessian was not computed,

   + Checking for singular values when the eigenvectors contained only
     one column did not work correctly,


 - Model comparison

   + Selection of the `modifiedharmonicmean` estimator was broken,


 - Optimal Simple Rules

   + When covariances were specified, variables that only entered with
     their variance and no covariance term obtained a wrong weight,
     resulting in wrong results,

   + Results reported for stochastic simulations after `osr` were based
     on the last parameter vector encountered during optimization,
     which does not necessarily coincide with the optimal parameter
     vector,

   + Using only one (co)variance in the objective function resulted in crashes,

   + For models with non-stationary variables the objective function was computed wrongly.


 - Ramsey policy

   + If a Lagrange multiplier appeared in the model with a lead or a lag
     of more than one period, the steady state could be wrong.

   + When using an external steady state file, incorrect steady states
     could be accepted,

   + When using an external steady state file with more than one
     instrument, Dynare crashed,

   + When using an external steady state file and running `stoch_simul`
     after `ramsey_planner`, an incorrect steady state was used,

   + When the number of instruments was not equal to the number of
     omitted equations, Dynare crashed with a cryptic message,

   + The `planner_objective` accepted `varexo`, but ignored them for computations,


 - Shock decomposition

   + Did not work with the `parameter_set=calibration` option if an
     `estimated_params` block is present,

   + Crashed after MLE.


 - Perfect foresight models

   + The perfect foresight solver could accept a complex solution
     instead of continuing to look for a real-valued one,

   + The `initval_file` command only accepted column and not row vectors,

   + The `initval_file` command did not work with Excel files,

   + Deterministic simulations with one boundary condition crashed in
     `solve_one_boundary` due to a missing underscore when passing
     `options_.simul.maxit`,

   + Deterministic simulation with exogenous variables lagged by more
     than one period crashed,

   + Termination criterion `maxit` was hard-coded for `solve_algo=0`
     and could no be changed,

   + When using `block`/`bytecode`, relational operators could not be enforced,

   + When using `block` some exceptions were not properly handled,
     leading to code crashes,

   + Using `periods=1` crashed the solver (bug only partially fixed).


 - Smoothing

   + The univariate Kalman smoother returned wrong results when used
     with correlated measurement error,

   + The diffuse smoother sometimes returned linear combinations of the
     smoothed stochastic trend estimates instead of the original trend
     estimates.

 - Perturbation reduced form

   + In contrast to what is stated in the manual, the results of the
     unconditional variance decomposition were only stored in
     `oo_.gamma_y(nar+2)`, not in `oo_.variance_decomposition`,

   + Dynare could crash when the steady state could not be computed
     when using the `loglinear` option,

   + Using `bytcode` when declared exogenous variables were not
     used in the model leaded to crashes in stochastic simulations,

   + Displaying decision rules involving lags of auxiliary variables of
     type 0 (leads>1) crashed.

   + The `relative_irf` option resulted in wrong output at `order>1` as
     it implicitly relies on linearity.


 - Displaying of the MH-history with the `internals` command crashed
   if parameter names did not have same length.

 - Dynare crashed when the user-defined steady state file returned an
   error code, but not an conformable-sized steady state vector.

 - Due to a bug in `mjdgges.mex` unstable parameter draws with
   eigenvalues up to 1+10⁻⁶ could be accepted as stable for the
   purpose of the Blanchard-Kahn conditions, even if `qz_criterium<1`.

 - The `use_dll` option on Octave for Windows required to pass a
   compiler flag at the command line, despite the manual stating this
   was not necessary.

 - Dynare crashed for models with `block` option if the Blanchard-Kahn
   conditions were not satisfied instead of generating an error
   message.

 - The `verbose` option did not work with `model(block)`.

 - When falsely specifying the `model(linear)` for nonlinear models,
   incorrect steady states were accepted instead of aborting.

 - The `STEADY_STATE` operator called on model local variables
   (so-called pound variables) did not work as expected.

 - The substring operator in macro-processor was broken. The
   characters of the substring could be mixed with random characters
   from the memory space.

 - Block decomposition could sometimes cause the preprocessor to crash.

 - A bug when external functions were used in model local variables
   that were contained in equations that required auxiliary
   variable/equations led to crashes of MATLAB.

 - Sampling from the prior distribution for an inverse gamma II
   distribution when `prior_trunc>0` could result in incorrect
   sampling.

 - Sampling from the prior distribution for a uniform distribution
   when `prior_trunc>0` was ignoring the prior truncation.

 - Conditional forecasts were wrong when the declaration of endogenous
   variables was not preceeding the declaration of the exogenous
   variables and parameters.



Announcement for Dynare 4.4.3 (on 2014-07-31)
=============================================

We are pleased to announce the release of Dynare 4.4.3.

This is a bugfix release.

The Windows packages are already available for download at:
<http://www.dynare.org/download/dynare-stable>.

The Mac and GNU/Linux packages (for Debian and Ubuntu LTS) should follow soon.

This release is compatible with MATLAB versions 7.3 (R2006b) to 8.2 (R2013b)
and with GNU Octave versions 3.6 to 3.8.

Here is a list of the problems identified in version 4.4.2 and that have been
fixed in version 4.4.3:

 - When loading a dataset in XLS, XLSX or CSV format, the first
   observation was discarded.

 - Reading data in an Excel-file with only one variable was leading
   to a crash.

 - When using the `k_order_perturbation` option (which is implicit at
   3rd order) without the `use_dll` option, crashes or unexpected
   behavior could happen if some 2nd or 3rd derivative evaluates to
   zero (while not being symbolically zero)

 - When using external function, Ramsey policy could crash or return
   wrong results.

 - For Ramsey policy, the equation numbers associated with the
   Lagrange multipliers stored in `M_.aux_vars` were erroneously one too
   low

 - When updating deep parameters in the steady state file, the changes
   were not fully taken into account (this was only affecting the
   Ramsey policy).

 - When using external functions and the bytecode option, wrong
   results were returned (if second order derivates of the external
   functions were needed).

 - The confidence level for computations in estimation, `conf_sig` could
   not be changed and was fixed at 0.9. The new option `mh_conf_sig` is
   now used to set this interval

 - Conditional forecasts with non-diagonal covariance matrix used an
   incorrect decomposition of the covariance matrix. A Cholesky
   factorization is used.

 - Option `geweke_interval` was not effective, Dynare always defaulted
   to the standard value.

 - The `mode_file` option lacked backward compatibility with older
   Dynare versions.

 - Loading an `mh_mode` file with the `mode_file` option was broken.

 - Using `identification` with `var_exo_det` leaded to crashes (the
   preprocessor now returns an error if they are used simultaneously)

 - The `identification` command did not print results if the initial
   parameter set was invalid and then crashed later on if the MC
   sample is bigger than 1

 - Inconsistencies between static and dynamic models leaded to crashes
   instead of error messages (only with block option).

 - The use of external functions crashed the preprocessor when the
   derivatives of the external function are explicitly called in the
   `model` block. The preprocessor now forbids the use of external
   functions derivates in the `model` block.

 - Using the block option when a variable does not appear in the
   current period crashed Dynare instead of providing an error
   message.


Announcement for Dynare 4.4.2 (on 2014-03-04)
=============================================

We are pleased to announce the release of Dynare 4.4.2.

This is a bugfix release.

The Windows packages are already available for download at:
<http://www.dynare.org/download/dynare-stable>.

The Mac and GNU/Linux packages (for Debian and Ubuntu LTS) should follow soon.

This release is compatible with MATLAB versions 7.3 (R2006b) to 8.2 (R2013b)
and with GNU Octave versions 3.6 to 3.8.

Here is a list of the problems identified in version 4.4.1 and that have been
fixed in version 4.4.2:

 - Geweke convergence diagnostics was computed on the wrong sample if `mh_drop`
   was not equal to the default of 0.5.

 - The `loglinear` option of `stoch_simul` was displaying the steady state of
   the original values, not the logged ones, and was producing incorrect
   simulations and simulated moments. Theoretical moments were unaffected.

 - The `optim` option of `estimation` (for setting options to `mode_compute`)
   was only working with at least MATLAB 8.1 (R2013a) or Octave 3.8.

 - For unit root models, theoretical HP filtered moments were sometimes
   erroneously displayed as NaN.

 - Specifying an endogenous variable twice after the `estimation` command would
   lead to a crash in the computation of moments.

 - Deterministic simulations were crashing on some models with more than one
   lead or one lag on exogenous variables.

 - Homotopy in stochastic extended path with order greater than 0 was not
   working correctly (during the homotopy steps the perfect foresight model
   solver was called instead of the stochastic perfect foresight model solver).

 - MCMC convergence diagnostics were not computed if `mh_replic` was less than
   2000; the test now relies on the total number of iterations (this only makes
   a difference if option `load_mh_file` is used).


Announcement for Dynare 4.4.1 (on 2014-01-17)
=============================================

We are pleased to announce the release of Dynare 4.4.1.

This release contains a few changes to the user interface and fixes various
bugs. It also adds compatibility with Octave 3.8.

The Windows packages are already available for download at:
<http://www.dynare.org/download/dynare-stable>.

The Mac and GNU/Linux packages (for Debian and Ubuntu) should follow soon.

All users are encouraged to upgrade.

This release is compatible with MATLAB versions 7.3 (R2006b) to 8.2 (R2013b) and
with GNU Octave versions 3.6 to 3.8.

* Changes to the user interface:

  - The syntax introduced in 4.4.0 for conditional forecast in a deterministic
    setup was removed, and replaced by a new one that is better suited to the
    task. More precisely, such deterministic forecasts are no longer done using
    the `conditional_forecast` command. The latter is replaced by a group of
    commands: `init_plan`, `basic_plan` and `flip_plan`. See the reference
    manual for more details.
 
  - Changes to the reporting module: option `annualAverages` to `addTable` has
    been removed (use option `tableDataRhs` to `addSeries` instead); option
    `vlineAfter` to `addTable` now also accepts a cell array.
 
  - Changes to the date and time series classes: implement broadcasting for
    operations (`+`,`-`,`*` and `/`) between `dseries` class and scalar or vectors; add
    the possibility of selecting an observation within a time series using a
    formatted string containing a date.
 
* Bugs and problems identified in version 4.4.0 and that have been fixed in
  version 4.4.1:

  - In MS-SBVAR, there was a bug preventing the computation of impulse responses
    on a constant regime.
 
  - Under Octave, after modifying the MOD file, the changes were not taken into
    account at the first Dynare run, but only at the second run.
 

Announcement for Dynare 4.4.0 (on 2013-12-16)
=============================================

We are pleased to announce the release of Dynare 4.4.0.

This major release adds new features and fixes various bugs.

The Windows packages are already available for download at:
<http://www.dynare.org/download/dynare-stable>.

The Mac and Debian/Ubuntu packages should follow soon.

All users are strongly encouraged to upgrade.

This release is compatible with MATLAB versions ranging from 7.3 (R2006b) to
8.2 (R2013b) and with GNU Octave version 3.6.

Here is the list of major user-visible changes:


* New major algorithms:

  - Extended path at order 1 and above, also known as “stochastic extended
    path”. This method is triggered by setting the `order` option of the
    `extended_path` command to a value greater than 0. Dynare will then use a
    Gaussian quadrature to take into account the effects of future uncertainty.
    The time series for the endogenous variables are generated by assuming that
    the agents believe that there will no more shocks after period t+order.
 
  - Alternative algorithms for computing decision rules of a stochastic model,
    based on the cycle reduction and logarithmic reduction algorithms. These
    methods are respectively triggered by giving `dr = cycle_reduction` or `dr
    = logarithmic_reduction` as an option to the `stoch_simul` command.
 
  - Pruning now works with 3rd order approximation, along the lines of
    *Andreasen, Fernández-Villaverde and Rubio-Ramirez (2013)*.
 
  - Computation of conditional forecast using an extended path method. This is
    triggered by the new option `simulation_type = deterministic` in the
    `conditional_forecast` command. In this case, the `expectation` command in
    the `conditional_forecast_paths` block has to be used to indicate the nature
    of expectations (whether shocks are a surprise or are perfectly
    anticipated).
 
  - Endogenous priors as in Christiano, Trabandt and Walentin (2011). Those are
    triggered by the new option `endogenous_prior` of the `estimation` command.
 

* Other algorithmic improvements:

  - New command `model_diagnostics` to perform various sanity checks on the
    model. Note: in the past, some users may have used a preliminary MATLAB
    function implementing this; the new command has the same syntax, except that
    you shouldn’t pass any argument to it.
 
  - Terminal conditions of perfect foresight simulations can now be specified in
    growth rates. More specifically, the new option `differentiate_forward_vars`
    of the `model` block will create auxiliary forward looking variables
    expressed in first differences or growth rates of the actual forward looking
    variables defined in the model. These new variables have obvious zero
    terminal conditions whatever the simulation context and this in many cases
    helps convergence of simulations.
 
  - Convergence diagnostics for single chain MCMC à la *Geweke (1992, 1999)*.
 
  - New optimizer for the posterior mode (triggered by `mode_compute=10`): it
    uses the simpsa algorithm, based on the combination of the non-linear
    simplex and simulated annealing algorithms and proposed by *Cardoso, Salcedo
    and Feyo de Azevedo (1996)*.
 
  - The automatic detrending engine has been extended to work on models written
    in logs. The corresponding trend variable type is `log_trend_var`, and the
    corresponding deflator type is `log_deflator`.
 

* New features in the user interface:

  - New set of functions for easily creating PDF reports including figures and
    tables. See the “Reporting” section in the reference manual for more
    details.
 
  - New MATLAB/Octave classes for handling time series. See the “Time series”
    section in the reference manual for more details.
 
  - Datafiles in CSV format can now be used for estimation.
 
  - New macro processor `length` operator, returns the length of an array.
 
  - New option `all_values_required` of `initval` and `endval` blocks: enforces
    initialization of all endogenous and exogenous variables within the block.
 
  - Option `ar` can now be given to the `estimation` command.
 
  - New options `nograph`, `nointeractive` and `nowarn` to the `dynare` command,
    for a better control of what is displayed.
 
  - New option `nostrict` to the `dynare` command, for allowing Dynare to
    continue processing when there are more endogenous variables than equations
    or when an undeclared symbol is assigned in `initval` or `endval`.
 
  - The information on MCMC acceptance rates, seeds, last log posterior
    likelihood, and last parameter draw are now saved on the disk and can
    be displayed with `internals --display-mh-history` or loaded into the
    workspace with `internals --load-mh-history`.
 
  - New options `mode_check_neighbourhood_size`, `mode_check_symmetric_plots`
    and `mode_check_number_of_points`, for a better control of the diagnostic
    plots.
 
  - New option `parallel_local_files` of `model` block, for transferring extra
    files during parallel computations.
 
  - New option `clock` of `set_dynare_seed`, for setting a different seed at
    each run.
 
  - New option `qz_zero_threshold` of the `check`, `stoch_simul` and
    `estimation` commands, for a better control of the situation where a
    generalized eigenvalue is close to 0/0.
 
  - New `verbatim` block for inclusion of text that should pass through the
    preprocessor and be placed as is in the `modfile.m` file.
 
  - New option `mcmc_jumping_covariance` of the `estimation` command, for a
    better control of the covariance matrix used for the proposal density of the
    MCMC sampler.
 
  - New option `use_calibration` of the `estimated_params_init`, for using the
    calibration of deep parameters and the elements of the covariance matrix
    specified in the `shocks` block as starting values for the estimation.
 
  - New option `save_draws` of the `ms_simulation` command.
 
  - New option `irf_plot_threshold` of the `stoch_simul` and `estimation`
    commands, for a better control of the display of IRFs which are almost nil.
 
  - New option `long_name` for endogenous, exogenous and parameter declarations,
    which can be used to declare a long name for variables. That long name can
    be programmatically retrieved in `M_.endo_names_long`.
 

* Miscellaneous changes

  - The deciles of some posterior moments were erroneously saved in a field
    `Distribution` under `oo_`. This field is now called `deciles`, for
    consistency with other posterior moments and with the manual. Similarly, the
    fields `Mean`, `Median`, `HPDsup`, `HPDinf`, and `Variance` are now
    consistently capitalized.
 
  - The console mode now implies the `nodisplay` option.
 

* Bugs and problems identified in version 4.3.3 and that have been fixed in
  version 4.4.0:

  - In an `endval` block, auxiliary variables were not given the right value.
    This would not result in wrong results, but could prevent convergence of
    the steady state computation.
 
  - Deterministic simulations with `stack_solve_algo=0` (the default value) were
    crashing if some exogenous had a lag strictly greater than 1.
 
  - When using the `mode_file` option, the initial estimation checks were not
    performed for the loaded mode, but for the original starting values. Thus,
    potential prior violations by the mode only appeared during estimation,
    leading to potentially cryptic crashes and error messages.
 
  - If a shock/measurement error variance was set to 0 in calibration, the
    correlation matrix featured a 0 instead of a 1 on the diagonal, leading to
    wrong estimation results.
 
  - In the presence of calibrated covariances, estimation did not enforce
    positive definiteness of the covariance matrix.
 
  - Estimation using the `diffuse_filter` option together with the univariate
    Kalman filter and a diagonal measurement error matrix was broken.
 
  - A purely backward model with `k_order_solver` was leading to crashes of
    MATLAB/Octave.
 
  - Non-linear estimation was not skipping the specified presample when
    computing the likelihood.
 
  - IRFs and theoretical moments at order > 1 were broken for purely
    forward-looking models.
 
  - Simulated moments with constant variables was leading to crashes when
    displaying autocorrelations.
 
  - The `osr` command was sometimes crashing with cryptic error messages because
    of some unaccounted error codes returned from a deeper routine.
 
  - The check for stochastic singularity during initial estimation checks was
    broken.
 
  - Recursive estimation starting with the pathological case of `nobs=1` was
    crashing.
 
  - Conditional variance decomposition within or after estimation was crashing
    when at least one shock had been calibrated to zero variance.
 
  - The `estimated_params_init` and `estimated_params_bounds` blocks were broken
    for correlations.
 
  - The `filter_step_ahead` option was not producing any output in Bayesian
    estimation.
 
  - Deterministic simulations were sometimes erroneously indicating convergence
    although the residuals were actually NaN or Inf.
 
  - Supplying a user function in the `mode_compute` option was leading to
    a crash.
 
  - Deterministic simulation of models without any exogenous variable was
    crashing.
 
  - The MS-SBVAR code was not updating files between runs on Windows. This means
    that if a MOD file was updated between runs in the same folder and a
    `file_tag` was not changed, then the results would not change.
 
  - The `ramsey_policy` command was not putting in `oo_.planner_objective_value`
    the value of the planner objective at the optimum.
 

* References:

  - Andreasen, Martin M., Jesús Fernández-Villaverde, and Juan Rubio-Ramirez
    (2013): “The Pruned State-Space System for Non-Linear DSGE Models: Theory
    and Empirical Applications,” *NBER Working Paper*, 18983
 
  - Cardoso, Margarida F., R. L. Salcedo and S. Feyo de Azevedo (1996): “The
    simplex simulated annealing approach to continuous non-linear optimization,”
    *Computers chem. Engng*, 20(9), 1065-1080
 
  - Christiano, Lawrence J., Mathias Trabandt and Karl Walentin (2011):
    “Introducing financial frictions and unemployment into a small open economy
    model,” *Journal of Economic Dynamics and Control*, 35(12), 1999-2041
 
  - Geweke, John (1992): “Evaluating the accuracy of sampling-based approaches
    to the calculation of posterior moments,” in J.O. Berger, J.M. Bernardo,
    A.P. Dawid, and A.F.M. Smith (eds.) *Proceedings of the Fourth Valencia
    International Meeting on Bayesian Statistics*, pp. 169-194, Oxford University
    Press
 
  - Geweke, John (1999): “Using simulation methods for Bayesian econometric
    models: Inference, development and communication,” *Econometric Reviews*,
    18(1), 1-73
 

Announcement for Dynare 4.3.3 (on 2013-04-12)
=============================================

We are pleased to announce the release of Dynare 4.3.3.

This is a bugfix release.

The Windows packages are already available for download at:
<http://www.dynare.org/download/dynare-stable>.

The Mac and GNU/Linux packages (for Debian and Ubuntu) should follow soon.

All users are encouraged to upgrade.

The new release is compatible with MATLAB versions ranging from 7.0 (R14) to
8.1 (R2013a) and with GNU Octave versions ranging from 3.2 to 3.6.

Here is a list of the problems identified in version 4.3.2 and that have been
fixed in version 4.3.3:

 - Estimation with measurement errors was wrong if a correlation between two
   measurement errors was calibrated

 - Option `use_dll` was broken under Windows

 - Degenerate case of purely static models (no leads/no lags) were not
   correctly handled

 - Deterministic simulations over a single period were not correctly done

 - The sensitivity call `dynare_sensitivity(identification=1,morris=2)` was
   buggy when there are no shocks estimated

 - Calls to `shock_decomposition` after using `selected_variables_only` option
   fail

 - Sometimes, only the last open graph was saved, leading to missing and
   duplicate EPS/PDF graphs

 - Forecasting after maximum likelihood estimation when not forecasting at
   least one observed variables (`var_obs`) was leading to crashes

 - Some functionalities were crashing with MATLAB 8.1/R2013a (bytecode,
   MS-SBVAR)

 - Sometimes only the first order autocorrelation of `moments_varendo` was
   saved instead of all up to the value of `ar` option


Announcement for Dynare 4.3.2 (on 2013-01-18)
=============================================

We are pleased to announce the release of Dynare 4.3.2.

This is a bugfix release.

The Windows packages are already available for download at:
<http://www.dynare.org/download/dynare-stable>.

The Mac and GNU/Linux packages (for Debian and Ubuntu) should follow soon.

All users are encouraged to upgrade.

The new release is compatible with MATLAB versions ranging from 7.0 (R14) to
8.0 (R2012b) and with GNU Octave versions ranging from 3.2 to 3.6.

Here is a list of the problems identified in version 4.3.1 and that have been
fixed in version 4.3.2:

 - Computation of posterior distribution of unconditional variance
   decomposition was sometimes crashing (only for very large models)

 - Estimation with `mode_compute=6` was sometimes crashing

 - Derivative of `erf()` function was incorrect

 - The `check` command was not setting `oo_.dr.eigval` unless `stoch_simul` was
   also used

 - Computation of conditional forecast when the constraint is only on
   one period was buggy

 - Estimation with `mode_compute=3` was crashing under Octave


Announcement for Dynare 4.3.1 (on 2012-10-10)
=============================================

We are pleased to announce the release of Dynare 4.3.1. This release adds a few
minor features and fixes various bugs.

The Windows and Mac packages are already available for download at:
<http://www.dynare.org/download/dynare-stable>.

The GNU/Linux packages (for Debian and Ubuntu) should follow soon.

All users are strongly encouraged to upgrade.

The new release is compatible with MATLAB versions ranging from 7.0 (R14) to
8.0 (R2012b) and with GNU Octave versions ranging from 3.2 to 3.6.

Here is the list of the main user-visible changes:


* New features in the user interface:

  - New `@#ifndef` directive in the macro-processor
 
  - Possibility of simultaneously specifying several output formats in the
    `graph_format` option
 
  - Support for XLSX files in `datafile` option of `estimation` and in
    `initval_file`
 

* Bugs and problems identified in version 4.3.0 and that have been fixed in
  version 4.3.1:

  - Shock decomposition was broken
 
  - The welfare computation with `ramsey_policy` was buggy when used in
    conjunction with `histval`
 
  - Estimation of models with both missing observations and measurement errors
    was buggy
 
  - The option `simul_replic` was broken
 
  - The macro-processor directive `@#ifdef` was broken
 
  - Identification with `max_dim_cova_group > 1` was broken for specially
    degenerate models (when parameter theta has pairwise collinearity of one
    with multiple other parameters, *i.e.* when all couples (θ,b), (θ,c),
    … (θ,d) have perfect collinearity in the Jacobian of the model)
 
  - The `parallel_test` option was broken
 
  - Estimation with correlated shocks was broken when the correlations were
    specified in terms of correlation and not in terms of co-variance
 
  - The Windows package was broken with MATLAB 7.1 and 7.2
 
  - When using `mode_compute=0` with a mode file generated using
    `mode_compute=6`, the value of option `mh_jscale` was not loaded
 
  - Using exogenous deterministic variables at 2nd order was causing a crash
 
  - The option `no_create_init` for the `ms_estimation` command was broken
 
  - Loading of datafiles with explicit filename extensions was not working
 
  - The preprocessor had a memory corruption problem which could randomly lead
    to crashes
 

Announcement for Dynare 4.3.0 (on 2012-06-15)
=============================================

We are pleased to announce the release of Dynare 4.3.0. This major release adds
new features and fixes various bugs.

The Windows and Mac packages are already available for download at:
<http://www.dynare.org/download/dynare-4.3>.

The GNU/Linux packages should follow soon.

All users are strongly encouraged to upgrade.

The new release is compatible with MATLAB versions ranging from 7.0 (R14) to
7.14 (R2012a) and with GNU Octave versions ranging from 3.2 to 3.6.

Here is the list of the main user-visible changes:


* New major algorithms:

  - Nonlinear estimation with a particle filter based on a second order
    approximation of the model, as in *Fernández-Villaverde and Rubio-Ramirez
    (2005)*; this is triggered by setting `order=2` in the `estimation` command
 
  - Extended path solution method as in *Fair and Taylor (1983)*; see the
    `extended_path` command
 
  - Support for Markov-Switching Structural Bayesian VARs (MS-SBVAR) along the
    lines of *Sims, Waggoner and Zha (2008)* (see the dedicated section in the
    reference manual)
 
  - Optimal policy under discretion along the lines of *Dennis (2007)*; see the
    `discretionary_policy` command
 
  - Identification analysis along the lines of *Iskrev (2010)*; see the
    `identification` command
 
  - The Global Sensitivity Analysis toolbox (*Ratto, 2008*) is now part of the
    official Dynare distribution
 

* Other algorithmic improvements:

  - Stochastic simulation and estimation can benefit from block decomposition
    (with the `block` option of `model`; only at 1st order)
 
  - Possibility of running smoother and filter on a calibrated model; see the
    `calib_smoother` command
 
  - Possibility of doing conditional forecast on a calibrated model; see the
    `parameter_set=calibration` option of the `conditional_forecast` command
 
  - The default algorithm for deterministic simulations has changed and is now
    based on sparse matrices; the historical algorithm (*Laffargue, Boucekkine
    and Juillard*) is still available under the `stack_solve_algo=6` option of the
    `simul` command
 
  - Possibility of using an analytic gradient for the estimation; see the
    `analytic_derivation` option of the `estimation` command
 
  - Implementation of the Nelder-Mead simplex based optimization routine for
    computing the posterior mode; available under the `mode_compute=8` option of
    the `estimation` command
 
  - Implementation of the CMA Evolution Strategy algorithm for computing the
    posterior mode; available under the `mode_compute=9` option of the
    `estimation` command
 
  - New solvers for Lyapunov equations which can accelerate the estimation of
    large models; see the `lyapunov` option of the `estimation` command
 
  - New solvers for Sylvester equations which can accelerate the resolution of
    large models with block decomposition; see the `sylvester` option of the
    `stoch_simul` and `estimation` commands
 
  - The `ramsey_policy` command now displays the planner objective value
    function under Ramsey policy and stores it in `oo_.planner_objective_value`
 
  - Theoretical autocovariances are now computed when the `block` option is
    present
 
  - The `linear` option is now compatible with the `block` and `bytecode`
    options
 
  - The `loglinear` option now works with purely backward or forward models at
    first order
 

* New features in the user interface:

  - New mathematical primitives allowed in model block: `abs()`, `sign()`
 
  - The behavior with respect to graphs has changed:
 
     + By default, Dynare now displays graphs and saves them to disk in EPS
       format only
 
     + The format can be changed to PDF or FIG with the new `graph_format`
       option
 
     + It is possible to save graphs to disk without displaying them with the
       new `nodisplay` option
 
  - New `nocheck` option to the `steady` command: tells not to check the steady
    state and accept values given by the user (useful for models with unit
    roots)
 
  - A series of deterministic shocks can be passed as a pre-defined vector in
    the `values` statement of a `shocks` block
 
  - New option `sub_draws` in the `estimation` command for controlling the
    number of draws used in computing the posterior distributions of various
    objects
 
  - New macroprocessor command `@#ifdef` for testing if a macro-variable is
    defined
 
  - New option `irf_shocks` of the `stoch_simul` command, to allow IRFs to be
    created only for certain exogenous variables
 
  - In the parallel engine, possibility of assigning different weights to nodes
    in the cluster and of creating clusters comprised of nodes with different
    operating systems (see the relevant section in the reference manual)
 
  - It is now possible to redefine a parameter in the `steady_state_model` block
    (use with caution)
 
  - New option `maxit` in the `simul` and `steady` commands to determine the
    maximum number of iterations of the nonlinear solver
 
  - New option `homotopy_force_continue` in the `steady` command to control the
    behavior when a homotopy fails
 
  - Possibility of globally altering the defaults of options by providing a file
    in the `GlobalInitFile` field of the configuration file (use with caution)
 
  - New option `nolog` to the `dynare` command line to avoid creating a logfile
 
  - New option `-D` to the `dynare` command line with for defining
    macro-variables
 

* Miscellaneous changes:

  - The `use_dll` option of `model` now creates a MEX file for the static model
    in addition to that for the dynamic model
 
  - The `unit_root_vars` command is now obsolete; use the `diffuse_filter`
    option of the `estimation` command instead
 
  - New option `--burn` to Dynare++ to discard initial simulation points
 
  - New top-level MATLAB/Octave command `internals` for internal documentation
    and unitary tests
 

* Bugs and problems identified in version 4.2.5 and that have been fixed in
  version 4.3.0:

  - Backward models with the `loglinear` option were incorrectly handled
 
  - Solving for hyperparameters of inverse gamma priors was sometimes crashing
 
  - The deterministic solver for purely forward models was broken
 
  - When running `estimation` or `identification` on models with non-diagonal
    structural error covariance matrices, while not simultaneously estimating
    the correlation between shocks (*i.e.* calibrating the correlation), the
    off-diagonal elements were incorrectly handled or crashes were occuring
 
  - When using the `prefilter` option, smoother plots were omitting the smoothed
    observables
 
  - In the rare case of entering and expression `x` as `x^(alpha-1)` with `x` being 0
    in steady state and alpha being a parameter equal to 2, the Jacobian was
    evaluating to 0 instead of 1
 
  - Setting the prior for shock correlations was failing if a lower bound was not
    explicitly specified
 

* References:

  - Dennis, Richard (2007): “Optimal Policy In Rational Expectations Models: New
    Solution Algorithms,” *Macroeconomic Dynamics*, 11(1), 31–55
 
  - Fair, Ray and John Taylor (1983): “Solution and Maximum Likelihood
    Estimation of Dynamic Nonlinear Rational Expectation Models,” *Econometrica*,
    51, 1169–1185
 
  - Fernández-Villaverde, Jesús and Juan Rubio-Ramirez (2005): “Estimating
    Dynamic Equilibrium Economies: Linear versus Nonlinear Likelihood,” *Journal
    of Applied Econometrics*, 20, 891–910
 
  - Iskrev, Nikolay (2010): “Local identification in DSGE models,” *Journal of
    Monetary Economics*, 57(2), 189–202
 
  - Ratto, Marco (2008): “Analysing DSGE models with global sensitivity
    analysis,” *Computational Economics*, 31, 115–139
 
  - Sims, Christopher A., Daniel F. Waggoner and Tao Zha (2008): “Methods for
    inference in large multiple-equation Markov-switching models,” *Journal of
    Econometrics*, 146, 255–274
 


Announcement for Dynare 4.2.5 (on 2012-03-14)
=============================================

We are pleased to announce the release of Dynare 4.2.5.

This is a bugfix release.

The Windows package for the new release is already available for download at
the official Dynare website <http://www.dynare.org>. The Mac and Linux packages
should follow soon.

All users are strongly encouraged to upgrade.

The new release is compatible with MATLAB versions ranging from 7.0 (R14) to
7.14 (R2012a) and with GNU Octave versions ranging from 3.0 to 3.6.

Note that GNU Octave users under Windows will have to upgrade to GNU Octave
version 3.6.1 (MinGW). The Octave installer can be downloaded at:
<http://www.dynare.org/octave/Octave3.6.1_gcc4.6.2_20120303-setup.exe>.

Here is a non-exhaustive list of the problems identified in version 4.2.4 and
that have been fixed in version 4.2.5:

 * The MATLAB optimization toolbox was sometimes not correctly detected even
   when installed

 * Using the inverse gamma distribution with extreme hyperparameter values
   could lead to a crash

 * Various issues in the accelerated deterministic solver with block
   decomposition

 * Various issues in the parallelization engine

 * Compatibility issues with the Global Sensitivity Analysis toolbox

 * The Dynare++ binary was broken in the Windows package because of a missing
   dynamic library


Announcement for Dynare 4.2.4 (on 2011-12-02)
=============================================

We are pleased to announce the release of Dynare 4.2.4.

This is a bugfix release. It comes only a few days after the previous release,
because version 4.2.3 was affected by a critical bug (see below).

The Windows package for the new release is already available for download at
the official [Dynare website](http://www.dynare.org). The Mac and Linux packages
should follow soon.

All users are strongly encouraged to upgrade, especially those who have
installed the buggy 4.2.3 release.

The new release is compatible with MATLAB versions ranging from 7.0 (R14) to
7.13 (R2011b) and with GNU Octave versions ranging from 3.0 to 3.4.

Here is the list of the problems identified in version 4.2.3 and that have been
fixed in version 4.2.4:

 * Second order approximation was broken for most models, giving incorrect
   results (this problem only affects version 4.2.3, not previous versions)

 * Bayesian priors with inverse gamma distribution and very small variances
   were giving incorrect results in some cases

 * The `model_diagnostics` command was broken


Announcement for Dynare 4.2.3 (on 2011-11-30)
=============================================

We are pleased to announce the release of Dynare 4.2.3.

This is a bugfix release.

The Windows package is already available for download at the official
[Dynare website](http://www.dynare.org). The Mac and Linux packages
should follow soon.

All users are strongly encouraged to upgrade.

This release is compatible with MATLAB versions ranging from 7.0 (R14)
to 7.13 (R2011b) and with GNU Octave versions ranging from 3.0 to 3.4.

Here is a non-exhaustive list of the problems identified in version 4.2.2 and
that have been fixed in version 4.2.3:

 * `steady_state_model` was broken for lags higher than 2

 * `simult_.m` was not working correctly with `order=3` if `k_order_solver` had
   not been explicitly specified

 * `stoch_simul` with `order=3` and without `periods` option was reporting
   dummy theoretical moments

 * Under Octave, option `solve_algo=0` was causing crashes in `check` and
   `stoch_simul`

 * Identification module was broken

 * The test for singularity in the model reporting eigenvalues close to 0/0 was
   sometimes reporting false positives

 * The `conditional_variance_decomposition` option was not working if one
   period index was 0. Now, Dynare reports an error if the periods are not
   strictly positive.

 * Second order approximation was buggy if one variable was not present at the
   current period


Announcement for Dynare 4.2.2 (on 2011-10-04)
=============================================

We are pleased to announce the release of Dynare 4.2.2.

This is a bugfix release.

The Windows package is already available for download at the official
[Dynare website](http://www.dynare.org). The Mac and Linux packages
should follow soon.

All users are strongly encouraged to upgrade.

This release is compatible with MATLAB versions ranging from 7.0 (R14)
to 7.13 (R2011b) and with GNU Octave versions ranging from 3.0 to 3.4.

Here is a list of the problems identified in version 4.2.1 and that have
been fixed in version 4.2.2:

 * The secondary rank test following the order test of the Blanchard and
   Kahn condition was faulty and almost never triggered

 * The variance prior for BVAR “à la Sims” with only one lag was
   inconsistent.  The solution implemented consists of adding one extra
   observation in the presample used to compute the prior; as a
   consequence, the numerical results for all estimations will be
   slightly different in future releases (thanks to Marek Jarociński for
   spotting this)

 * The `conditional_forecast` command was buggy: it was always using the
   posterior mode, whatever the value of the `parameter_set` option

 * `STEADY_STATE` was not working correctly with certain types of
   expressions (the priority of the addition and substraction operators
   was incorrectly handled)

 * With the `block` option of `model`, the preprocessor was failing on
   expressions of the form `a^b` (with no endogenous in `a` but an
   endogenous in `b`)

 * Some native MATLAB statements were not correctly passed on to MATLAB
   (*e.g.* `x = { 'foo' 'bar' }` )

 * `external_function` was crashing in some circumstances

 * The lambda parameter for HP filter was restricted to integer values
   for no good reason

 * The `load_mh_file` option of `estimation` was crashing under Octave
   for Windows (MinGW version)

 * Computation of steady state was failing on model contains auxiliary
   variables created by leads or lags larger than 2 or by of the
   `EXPECTATION` operator

 * Compilation of MEX files for MATLAB was failing with GCC 4.6


Announcement for Dynare 4.2.1 (on 2011-05-24)
=============================================

We are pleased to announce the release of Dynare 4.2.1.

Many bugs have been fixed since the previous release. The reference
manual has also been improved: new contents has been added at various
places, the structure has been improved, an index of functions and
variables has been added, the PDF/HTML rendering has been improved.

The Windows package is already available for download at the official
[Dynare website](http://www.dynare.org). The Mac and Linux packages should follow soon.

All users are strongly encouraged to upgrade.

This release is compatible with MATLAB versions ranging from 7.0 (R14)
to 7.12 (R2011a) and with GNU Octave versions ranging from 3.0 to 3.4.

Here is a list of the main bugfixes since version 4.2.0:

 * The `STEADY_STATE` operator has been fixed

 * Problems with MATLAB 7.3 (R2006b) and older have been fixed

 * The `partial_information` option of `stoch_simul` has been fixed

 * Option `conditional_variance_decomposition` of `stoch_simul` and
   `estimation` has been fixed

 * Automatic detrending now works in conjunction with the `EXPECTATION`
   operator

 * Percentage signs inside strings in MATLAB statements (like `disp('%
   This is not a comment %')`) now work

 * Beta prior with a very small standard deviation now work even if you
   do not have the MATLAB Statistical toolbox

 * External functions can now been used in assignment of model local
   variables

 * `identification` command has been fixed

 * Option `cova_compute` of `estimation` command has been fixed

 * Random crashes with 3rd order approximation without `use_dll` option
   have been eliminated


Announcement for Dynare 4.2.0 (on 2011-02-15)
=============================================

We are pleased to announce the release of Dynare 4.2.0.

This major release adds new features and fixes various bugs.

The Windows package is already available for download. The Mac and Linux
packages should follow soon.

All users are strongly encouraged to upgrade.

This release is compatible with MATLAB versions ranging from 6.5 (R13) to 7.11
(R2010b) and with GNU Octave versions 3.0.x and 3.2.x (support for GNU Octave
3.4.x is not complete and will be added in the next minor release).

Here is the list of major user-visible changes:

* New solution algorithms:

  - Pruning for second order simulations has been added, as described in *Kim,
    Kim, Schaumburg and Sims (2008)*. It is triggered by option `pruning` of
    `stoch_simul` (only 2nd order, not available at 3rd order).

  - Models under partial information can be solved, as in *Pearlman, Currie and
    Levine (1986)*. See <http://www.dynare.org/DynareWiki/PartialInformation>.

  - New nonlinear solvers for faster deterministic simulations and steady state
    computation. See
    <http://www.dynare.org/DynareWiki/FastDeterministicSimulationAndSteadyStateComputation>.

* Dynare can now use the power of multi-core computers or of a cluster of
  computer using parallelization. See
  <http://www.dynare.org/DynareWiki/ParallelDynare>.

* New features in the user interface:

  - A steady state file can now be automatically generated, provided that the
    model can be solved analytically, and that the steady state as a function
    of the parameters is declared with the new `steady_state_model` command.
    See the entry for `steady_state_model` in the reference manual for more
    details and an example.

  - For non-stationary models, Dynare is now able of automatically removing
    trends in all the equations: the user writes the equations in
    non-stationary form and declares the deflator of each variable. Then Dynare
    perform a check to determine if the proposed deflators are compatible with
    balanced growth path, and, if yes, then it computes the detrended
    equations. See <http://www.dynare.org/DynareWiki/RemovingTrends>.

  - It is now possible to use arbitrary functions in the model block. See
    <http://www.dynare.org/DynareWiki/ExternalFunctions>.

* Other minor changes to the user interface:

  - New primitives allowed in model block: `normpdf()`, `erf()`

  - New syntax for DSGE-VAR. See <http://www.dynare.org/DynareWiki/DsgeVar>.

  - Syntax of deterministic shocks has changed: after the values keyword,
    arbitrary expressions must be enclosed within parentheses (but numeric
    constants are still accepted as is)

* Various improvements:

  - Third order simulations now work without the `use_dll` option:
    installing a C++ compiler is no longer necessary for 3rd order

  - The HP filter works for empirical moments (previously it was only available
    for theoretical moments)

  - `ramsey_policy` now displays the planner objective value function under
    Ramsey policy and stores it in `oo_.planner_objective_value`

  - Estimation: if the `selected_variables_only` option is present, then the
    smoother will only be run on variables listed just after the estimation
    command

  - Estimation: in the `shocks` block, it is now possible to calibrate
    measurement errors on endogenous variables (using the same keywords than
    for calibrating variance/covariance matrix of exogenous shocks)

  - It is possibile to choose the parameter set for shock decomposition. See
    <http://www.dynare.org/DynareWiki/ShockDecomposition>.

  - The diffuse filter now works under Octave

  - New option `console` on the Dynare command-line: use it when running Dynare
    from the console, it will replace graphical waitbars by text waitbars for
    long computations

  - Steady option `solve_algo=0` (uses `fsolve()`) now works under Octave

* For Emacs users:

  - New Dynare mode for Emacs editor (contributed by Yannick Kalantzis)

  - Reference manual now available in Info format (distributed with
    Debian/Ubuntu packages)

* Miscellaneous:

  - Deterministic models: leads and lags of two or more on endogenous
    variables are now substituted by auxiliary variables; exogenous variables
    are left as is. See <http://www.dynare.org/DynareWiki/AuxiliaryVariables>.

* References:

  - Kim, J., S. Kim, E. Schaumburg and C.A. Sims (2008), “Calculating and using
    second-order accurate solutions of discrete time dynamic equilibrium
    models,” *Journal of Economic Dynamics and Control*, 32(11), 3397-3414

  - Pearlman J., D. Currie and P. Levine (1986), “Rational expectations models
    with partial information,” *Economic Modelling*, 3(2), 90-105

