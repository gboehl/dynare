function [dataset_, dataset_info, xparam1, hh, M_, options_, oo_, estim_params_,bayestopt_, bounds] = dynare_estimation_init(var_list_, dname, gsa_flag, M_, options_, oo_, estim_params_, bayestopt_)

% Performs initialization tasks before estimation or global sensitivity analysis
%
% INPUTS
%   var_list_:      selected endogenous variables vector
%   dname:          alternative directory name
%   gsa_flag:       flag for GSA operation (optional)
%   M_:             structure storing the model information
%   options_:       structure storing the options
%   oo_:            structure storing the results
%   estim_params_:  structure storing information about estimated
%                   parameters
%   bayestopt_:     structure storing information about priors

% OUTPUTS
%   dataset_:       the dataset after required transformation
%   dataset_info:   Various informations about the dataset (descriptive statistics and missing observations).
%   xparam1:        initial value of estimated parameters as returned by
%                   set_prior() or loaded from mode-file
%   hh:             hessian matrix at the loaded mode (or empty matrix)
%   M_:             structure storing the model information
%   options_:       structure storing the options
%   oo_:            structure storing the results
%   estim_params_:  structure storing information about estimated
%                   parameters
%   bayestopt_:     structure storing information about priors
%   bounds:         structure containing prior bounds
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2003-2023 Dynare Team
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

hh = [];
xparam1 = [];

if isempty(gsa_flag)
    gsa_flag = false;
else
    % Decide if a DSGE or DSGE-VAR has to be estimated.
    if ~isempty(strmatch('dsge_prior_weight', M_.param_names))
        options_.dsge_var = 1;
    end
    if isempty(var_list_)
        var_list_ = check_list_of_variables(options_, M_, var_list_);
        options_.varlist = var_list_;
    end
    if gsa_flag
        % Get the list of the endogenous variables for which posterior statistics wil be computed.
        options_.varlist = var_list_;
    else
        % This was done in dynare_estimation_1
    end
end

if options_.dsge_var && options_.presample~=0
    error('DSGE-VAR does not support the presample option.')
end

% Test if observed variables are declared.
if ~isfield(options_,'varobs')
    error('VAROBS statement is missing!')
end

% Checks on VAROBS
check_varobs_are_endo_and_declared_once(options_.varobs,M_.endo_names);
% Set the number of observed variables.
options_.number_of_observed_variables = length(options_.varobs);

if options_.discretionary_policy
    if options_.order>1
        error('discretionary_policy does not support order>1');
    else
        M_=discretionary_policy_initialization(M_,options_);
    end
end

% Check the perturbation order for pruning (k order perturbation based nonlinear filters are not yet implemented for k>3).
if options_.order>3 && options_.particle.pruning
    error('Higher order nonlinear filters are not compatible with pruning option.')
end

% analytical derivation is not yet available for kalman_filter_fast
if options_.analytic_derivation && options_.fast_kalman_filter
    error(['estimation option conflict: analytic_derivation isn''t available ' ...
           'for fast_kalman_filter'])
end

% fast kalman filter is only available with kalman_algo == 0,1,3
if options_.fast_kalman_filter && ~ismember(options_.kalman_algo, [0,1,3])
    error(['estimation option conflict: fast_kalman_filter is only available ' ...
               'with kalman_algo = 0, 1 or 3'])
end

% Set options_.lik_init equal to 3 if diffuse filter is used or kalman_algo refers to a diffuse filter algorithm.
if isequal(options_.diffuse_filter,1) || (options_.kalman_algo>2)
    if isequal(options_.lik_init,2)
        error('options diffuse_filter, lik_init and/or kalman_algo have contradictory settings')
    else
        options_.lik_init = 3;
    end
end

options_=select_qz_criterium_value(options_);

% Set options related to filtered variables.
if  isequal(options_.filtered_vars,0) && ~isempty(options_.filter_step_ahead)
    options_.filtered_vars = 1;
end
if ~isequal(options_.filtered_vars,0) && isempty(options_.filter_step_ahead)
    options_.filter_step_ahead = 1;
end
if ~isequal(options_.filtered_vars,0) && isequal(options_.filter_step_ahead,0)
    options_.filter_step_ahead = 1;
end
if ~isequal(options_.filter_step_ahead,0)
    options_.nk = max(options_.filter_step_ahead);
end

% Set the name of the directory where (intermediary) results will be saved.
if isempty(dname)
    M_.dname = M_.fname;
else
    M_.dname = dname;
end

% Set priors over the estimated parameters.
if ~isempty(estim_params_) && ~(isfield(estim_params_,'nvx') && (size(estim_params_.var_exo,1)+size(estim_params_.var_endo,1)+size(estim_params_.corrx,1)+size(estim_params_.corrn,1)+size(estim_params_.param_vals,1))==0)
    [xparam1,estim_params_,bayestopt_,lb,ub,M_] = set_prior(estim_params_,M_,options_);
end

if ~isempty(bayestopt_) && any(bayestopt_.pshape==0) && any(bayestopt_.pshape~=0)
    error('Estimation must be either fully ML or fully Bayesian. Maybe you forgot to specify a prior distribution.')
end
% Check if a _prior_restrictions.m file exists
if exist([M_.fname '_prior_restrictions.m'])
    options_.prior_restrictions.status = 1;
    options_.prior_restrictions.routine = str2func([M_.fname '_prior_restrictions']);
end

% Check that the provided mode_file is compatible with the current estimation settings.
if ~isempty(estim_params_) && ~(isfield(estim_params_,'nvx') && sum(estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.np)==0) && ~isempty(options_.mode_file) && ~options_.mh_posterior_mode_estimation
    [xparam1, hh] = check_mode_file(xparam1, hh, options_, bayestopt_);
end

%check for calibrated covariances before updating parameters
if ~isempty(estim_params_) && ~(isfield(estim_params_,'nvx') && sum(estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.np)==0)
    estim_params_=check_for_calibrated_covariances(estim_params_,M_);
end

%%read out calibration that was set in mod-file and can be used for initialization
xparam1_calib=get_all_parameters(estim_params_,M_); %get calibrated parameters
if ~any(isnan(xparam1_calib)) %all estimated parameters are calibrated
    estim_params_.full_calibration_detected=1;
else
    estim_params_.full_calibration_detected=0;
end
if options_.use_calibration_initialization %set calibration as starting values
    if ~isempty(bayestopt_) && all(bayestopt_.pshape==0) && any(all(isnan([xparam1_calib xparam1]),2))
        error('Estimation: When using the use_calibration option with ML, the parameters must be properly initialized.')
    else
        [xparam1,estim_params_]=do_parameter_initialization(estim_params_,xparam1_calib,xparam1); %get explicitly initialized parameters that have precedence to calibrated values
    end
end

if ~isempty(bayestopt_) && all(bayestopt_.pshape==0) && any(isnan(xparam1))
    error('ML estimation requires all estimated parameters to be initialized, either in an estimated_params or estimated_params_init-block ')
end

if ~isempty(estim_params_) && ~(all(strcmp(fieldnames(estim_params_),'full_calibration_detected'))  || (isfield(estim_params_,'nvx') && sum(estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.np)==0))
    if ~isempty(bayestopt_) && any(bayestopt_.pshape > 0)
        % Plot prior densities.
        if ~options_.nograph && options_.plot_priors
            plot_priors(bayestopt_,M_,estim_params_,options_)
        end
        % Set prior bounds
        bounds = prior_bounds(bayestopt_, options_.prior_trunc);
        bounds.lb = max(bounds.lb,lb);
        bounds.ub = min(bounds.ub,ub);
    else  % estimated parameters but no declared priors
          % No priors are declared so Dynare will estimate the model by
          % maximum likelihood with inequality constraints for the parameters.
        options_.mh_replic = 0;% No metropolis.
        bounds.lb = lb;
        bounds.ub = ub;
    end
    % Test if initial values of the estimated parameters are all between the prior lower and upper bounds.
    if options_.use_calibration_initialization
        try
            check_prior_bounds(xparam1,bounds,M_,estim_params_,options_,bayestopt_)
        catch
            e = lasterror();
            fprintf('Cannot use parameter values from calibration as they violate the prior bounds.')
            rethrow(e);
        end
    else
        check_prior_bounds(xparam1,bounds,M_,estim_params_,options_,bayestopt_)
    end
end

if isempty(estim_params_) || all(strcmp(fieldnames(estim_params_),'full_calibration_detected')) || (isfield(estim_params_,'nvx') && sum(estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.np)==0) % If estim_params_ is empty (e.g. when running the smoother on a calibrated model)
    if ~options_.smoother
        error('Estimation: the ''estimated_params'' block is mandatory (unless you are running a smoother)')
    end
    xparam1 = [];
    bayestopt_.jscale = [];
    bayestopt_.pshape = [];
    bayestopt_.name =[];
    bayestopt_.p1 = [];
    bayestopt_.p2 = [];
    bayestopt_.p3 = [];
    bayestopt_.p4 = [];
    bayestopt_.p5 = [];
    bayestopt_.p6 = [];
    bayestopt_.p7 = [];
    estim_params_.var_exo=[];
    estim_params_.var_endo=[];
    estim_params_.corrx=[];
    estim_params_.corrn=[];
    estim_params_.param_vals=[];
    estim_params_.nvx = 0;
    estim_params_.nvn = 0;
    estim_params_.ncx = 0;
    estim_params_.ncn = 0;
    estim_params_.np = 0;
    bounds.lb = [];
    bounds.ub = [];
end

% storing prior parameters in results
oo_.prior.mean = bayestopt_.p1;
oo_.prior.mode = bayestopt_.p5;
oo_.prior.variance = diag(bayestopt_.p2.^2);
oo_.prior.hyperparameters.first = bayestopt_.p6;
oo_.prior.hyperparameters.second = bayestopt_.p7;

% Is there a linear trend in the measurement equation?
if ~isfield(options_,'trend_coeffs') % No!
    bayestopt_.with_trend = 0;
else% Yes!
    bayestopt_.with_trend = 1;
    bayestopt_.trend_coeff = {};
    for i=1:options_.number_of_observed_variables
        if i > length(options_.trend_coeffs)
            bayestopt_.trend_coeff{i} = '0';
        else
            bayestopt_.trend_coeff{i} = options_.trend_coeffs{i};
        end
    end
end

% Get informations about the variables of the model.
dr = set_state_space(oo_.dr,M_);
oo_.dr = dr;
nstatic = M_.nstatic;          % Number of static variables.
npred = M_.nspred;             % Number of predetermined variables.
nspred = M_.nspred;            % Number of predetermined variables in the state equation.

%% Setting restricted state space (observed + predetermined variables)
% oo_.dr.restrict_var_list: location of union of observed and state variables in decision rules (decision rule order)
% bayestopt_.mfys: position of observables in oo_.dr.ys (declaration order)
% bayestopt_.mf0: position of state variables in restricted state vector (oo_.dr.restrict_var_list)
% bayestopt_.mf1: positions of observed variables in restricted state vector (oo_.dr.restrict_var_list order)
% bayestopt_.mf2: positions of observed variables in decision rules/expanded state vector (decision rule order)
% bayestopt_.smoother_var_list: positions of observed variables and requested smoothed variables in decision rules (decision rule order)
% bayestopt_.smoother_saved_var_list: positions of requested smoothed variables in bayestopt_.smoother_var_list
% bayestopt_.smoother_restrict_columns: positions of states in observed variables and requested smoothed variables in decision rules (decision rule order)
% bayestopt_.smoother_mf: positions of observed variables and requested smoothed variables in bayestopt_.smoother_var_list
var_obs_index_dr = [];
k1 = [];
for i=1:options_.number_of_observed_variables
    var_obs_index_dr = [var_obs_index_dr; strmatch(options_.varobs{i}, M_.endo_names(dr.order_var), 'exact')];
    k1 = [k1; strmatch(options_.varobs{i}, M_.endo_names, 'exact')];
end

k3 = [];
k3p = [];
if options_.selected_variables_only
    if options_.forecast > 0 && options_.mh_replic == 0 && ~options_.load_mh_file
        fprintf('\nEstimation: The selected_variables_only option is incompatible with classical forecasts. It will be ignored.\n')
        k3 = (1:M_.endo_nbr)';
        k3p = (1:M_.endo_nbr)';
    else
        for i=1:length(var_list_)
            k3 = [k3; strmatch(var_list_{i}, M_.endo_names(dr.order_var), 'exact')];
            k3p = [k3; strmatch(var_list_{i}, M_.endo_names, 'exact')];
        end
    end
else
    k3 = (1:M_.endo_nbr)';
    k3p = (1:M_.endo_nbr)';
end

% Define union of observed and state variables
k2 = union(var_obs_index_dr,[M_.nstatic+1:M_.nstatic+M_.nspred]', 'rows');
% Set restrict_state to postion of observed + state variables in expanded state vector.
oo_.dr.restrict_var_list = k2;
% set mf0 to positions of state variables in restricted state vector for likelihood computation.
[~,bayestopt_.mf0] = ismember([M_.nstatic+1:M_.nstatic+M_.nspred]',k2);
% Set mf1 to positions of observed variables in restricted state vector for likelihood computation.
[~,bayestopt_.mf1] = ismember(var_obs_index_dr,k2);
% Set mf2 to positions of observed variables in expanded state vector for filtering and smoothing.
bayestopt_.mf2  = var_obs_index_dr;
bayestopt_.mfys = k1;
[~,ic] = intersect(k2,nstatic+(1:npred)');
oo_.dr.restrict_columns = [ic; length(k2)+(1:nspred-npred)'];
bayestopt_.smoother_var_list = union(k2,k3);
[~,~,bayestopt_.smoother_saved_var_list] = intersect(k3,bayestopt_.smoother_var_list(:));
[~,ic] = intersect(bayestopt_.smoother_var_list,nstatic+(1:npred)');
bayestopt_.smoother_restrict_columns = ic;
[~,bayestopt_.smoother_mf] = ismember(var_obs_index_dr, bayestopt_.smoother_var_list);

if options_.analytic_derivation
    if options_.lik_init == 3
        error('analytic derivation is incompatible with diffuse filter')
    end
    options_.analytic_derivation = 1;
    if estim_params_.np || isfield(options_,'identification_check_endogenous_params_with_no_prior')
        % check if steady state changes param values
        M_local=M_;
        if isfield(options_,'identification_check_endogenous_params_with_no_prior')
            M_local.params = M_local.params*1.01; %vary parameters
        else
            M_local.params(estim_params_.param_vals(:,1)) = xparam1(estim_params_.nvx+estim_params_.ncx+estim_params_.nvn+estim_params_.ncn+1:end); %set parameters
            M_local.params(estim_params_.param_vals(:,1)) = M_local.params(estim_params_.param_vals(:,1))*1.01; %vary parameters
        end
        if options_.diffuse_filter || options_.steadystate.nocheck
            steadystate_check_flag = 0;
        else
            steadystate_check_flag = 1;
        end
        [~, params] = evaluate_steady_state(oo_.steady_state,[oo_.exo_steady_state; oo_.exo_det_steady_state],M_local,options_,steadystate_check_flag);
        change_flag=any(find(params-M_local.params));
        if change_flag
            skipline()
            if any(isnan(params))
                disp('After computing the steadystate, the following parameters are still NaN: '),
                disp(char(M_local.param_names(isnan(params))))
            end
            if any(find(params(~isnan(params))-M_local.params(~isnan(params))))
                disp('The steadystate file changed the values for the following parameters: '),
                disp(char(M_local.param_names(find(params(~isnan(params))-M_local.params(~isnan(params))))))
            end
            disp('The derivatives of jacobian and steady-state will be computed numerically'),
            disp('(re-set options_.analytic_derivation_mode= -2)'),
            options_.analytic_derivation_mode= -2;
        end
    end
end

% If jscale isn't specified for an estimated parameter, use global option options_.jscale, set to optimal value for a normal distribution by default.
% Note that check_posterior_sampler_options and mode_compute=6 may overwrite the setting
if isempty(options_.mh_jscale)
    options_.mh_jscale=2.38/sqrt(length(xparam1));
end
k = find(isnan(bayestopt_.jscale));
bayestopt_.jscale(k) = options_.mh_jscale;

% Build the dataset
if ~isempty(options_.datafile)
    [~,name] = fileparts(options_.datafile);
    if strcmp(name,M_.fname)
        error('Data-file and mod-file are not allowed to have the same name. Please change the name of the data file.')
    end
end

if isnan(options_.first_obs)
    options_.first_obs=1;
end
[dataset_, dataset_info] = makedataset(options_, options_.dsge_var*options_.dsge_varlag, gsa_flag);

%set options for old interface from the ones for new interface
if ~isempty(dataset_)
    options_.nobs = dataset_.nobs;
    if options_.endogenous_prior
        if ~isnan(dataset_info.missing.number_of_observations) && ~(dataset_info.missing.number_of_observations==0) %missing observations present
            if dataset_info.missing.no_more_missing_observations<dataset_.nobs-10
                fprintf('\ndynare_estimation_init: There are missing observations in the data.\n')
                fprintf('dynare_estimation_init: I am computing the moments for the endogenous prior only\n')
                fprintf('dynare_estimation_init: on the observations after the last missing one, i.e. %u.\n',dataset_info.missing.no_more_missing_observations)
            else
                fprintf('\ndynare_estimation_init: There are too many missing observations in the data.\n')
                fprintf('dynare_estimation_init: The endogenous_prior-option needs a consistent sample of \n')
                fprintf('dynare_estimation_init: at least 10 full observations at the end.\n')
                error('The endogenous_prior-option does not support your missing data.')
            end
        end
    end
end

% setting steadystate_check_flag option
if options_.diffuse_filter || options_.steadystate.nocheck
    steadystate_check_flag = 0;
else
    steadystate_check_flag = 1;
end

%check steady state at initial parameters
M_local = M_;
nvx = estim_params_.nvx;
ncx = estim_params_.ncx;
nvn = estim_params_.nvn;
ncn = estim_params_.ncn;
if estim_params_.np
    M_local.params(estim_params_.param_vals(:,1)) = xparam1(nvx+ncx+nvn+ncn+1:end);
end
[oo_.steady_state, params,info] = evaluate_steady_state(oo_.steady_state,[oo_.exo_steady_state; oo_.exo_det_steady_state],M_local,options_,steadystate_check_flag);

if info(1)
    fprintf('\ndynare_estimation_init:: The steady state at the initial parameters cannot be computed.\n')
    if options_.debug
        M_local.params=params;
        plist = list_of_parameters_calibrated_as_NaN(M_local);
        if ~isempty(plist)
            message = 'dynare_estimation_init:: Some of the parameters are NaN (' ;
            for i=1:length(plist)
                if i<length(plist)
                    message = [message, plist{i} ', '];
                else
                    message = [message, plist{i} ')'];
                end
            end
        end
        fprintf('%s\n',message)
        plist = list_of_parameters_calibrated_as_Inf(M_local);
        if ~isempty(plist)
            message = 'dynare_estimation_init:: Some of the parameters are Inf (' ;
            for i=1:length(plist)
                if i<size(plist)
                    message = [message, plist{i} ', '];
                else
                    message = [message, plist{i} ')'];
                end
            end
        end
        fprintf('%s\n',message)
    end
    print_info(info, 0, options_);
end

% If the steady state of the observed variables is non zero, set noconstant equal 0 ()
if (~options_.loglinear && all(abs(oo_.steady_state(bayestopt_.mfys))<1e-9)) || (options_.loglinear && all(abs(log(oo_.steady_state(bayestopt_.mfys)))<1e-9))
    options_.noconstant = 0; %identifying the constant based on just the initial parameter value is not feasible
else
    options_.noconstant = 0;
    % If the data are prefiltered then there must not be constants in the
    % measurement equation of the DSGE model or in the DSGE-VAR model.
    if options_.prefilter
        skipline()
        disp('You have specified the option "prefilter" to demean your data but the')
        disp('steady state of of the observed variables is non zero.')
        disp('Either change the measurement equations, by centering the observed')
        disp('variables in the model block, or drop the prefiltering.')
        error('The option "prefilter" is inconsistent with the non-zero mean measurement equations in the model.')
    end
end

%% get the non-zero rows and columns of Sigma_e and H

estim_params_= get_matrix_entries_for_psd_check(M_,estim_params_);

if options_.load_results_after_load_mh
    if ~exist([M_.dname filesep 'Output' filesep M_.fname '_results.mat'],'file')
        fprintf('\ndynare_estimation_init:: You specified the load_results_after_load_mh, but no _results.mat-file\n')
        fprintf('dynare_estimation_init:: was found. Results will be recomputed.\n')
        options_.load_results_after_load_mh=0;
    end
end

if options_.mh_replic || options_.load_mh_file
    [current_options, options_, bayestopt_] = check_posterior_sampler_options([], M_.fname, M_.dname, options_, bounds, bayestopt_);
    options_.posterior_sampler_options.current_options = current_options;
end

%% set heteroskedastic shocks
if options_.heteroskedastic_filter
    if options_.fast_kalman_filter
        error('estimation option conflict: "heteroskedastic_filter" incompatible with "fast_kalman_filter"')
    end
    if options_.analytic_derivation
        error(['estimation option conflict: analytic_derivation isn''t available ' ...
            'for heteroskedastic_filter'])
    end

    M_.heteroskedastic_shocks.Qvalue = NaN(M_.exo_nbr,options_.nobs+1);
    M_.heteroskedastic_shocks.Qscale = NaN(M_.exo_nbr,options_.nobs+1);

    for k=1:length(M_.heteroskedastic_shocks.Qvalue_orig)
        v = M_.heteroskedastic_shocks.Qvalue_orig(k);
        temp_periods=v.periods(v.periods<options_.nobs+options_.first_obs+1);
        temp_periods=temp_periods(temp_periods>=options_.first_obs);
        M_.heteroskedastic_shocks.Qvalue(v.exo_id, temp_periods-(options_.first_obs-1)) = v.value^2;
    end
    for k=1:length(M_.heteroskedastic_shocks.Qscale_orig)
        v = M_.heteroskedastic_shocks.Qscale_orig(k);
        temp_periods=v.periods(v.periods<options_.nobs+options_.first_obs+1);
        temp_periods=temp_periods(temp_periods>=options_.first_obs);
        M_.heteroskedastic_shocks.Qscale(v.exo_id, temp_periods-(options_.first_obs-1)) = v.scale^2;
    end

    if any(any(~isnan(M_.heteroskedastic_shocks.Qvalue) & ~isnan(M_.heteroskedastic_shocks.Qscale)))
        fprintf('\ndynare_estimation_init: With the option "heteroskedastic_shocks" you cannot define\n')
        fprintf('dynare_estimation_init: the scale and the value for the same shock \n')
        fprintf('dynare_estimation_init: in the same period!\n')
        error('Scale and value defined for the same shock in the same period with "heteroskedastic_shocks".')
    end
end

if (options_.occbin.likelihood.status && options_.occbin.likelihood.inversion_filter) || (options_.occbin.smoother.status && options_.occbin.smoother.inversion_filter)
    if isempty(options_.occbin.likelihood.IVF_shock_observable_mapping)
        options_.occbin.likelihood.IVF_shock_observable_mapping=find(diag(M_local.Sigma_e)~=0);
    else
        zero_var_shocks=find(diag(M_local.Sigma_e)==0);
        if any(ismember(options_.occbin.likelihood.IVF_shock_observable_mapping,zero_var_shocks))
            error('IVF-filter: an observable is mapped to a zero variance shock.')
        end
    end
end

if options_.occbin.smoother.status && options_.occbin.smoother.inversion_filter
    if ~isempty(options_.nk)
        fprintf('dynare_estimation_init: the inversion filter does not support filter_step_ahead. Disabling the option.\n')
        options_.nk=[];
    end
    if options_.filter_covariance
        fprintf('dynare_estimation_init: the inversion filter does not support filter_covariance. Disabling the option.\n')
        options_.filter_covariance=false;
    end
    if options_.smoothed_state_uncertainty
        fprintf('dynare_estimation_init: the inversion filter does not support smoothed_state_uncertainty. Disabling the option.\n')
        options_.smoothed_state_uncertainty=false;
    end
end

if options_.occbin.smoother.status || options_.occbin.likelihood.status
    if isfield(M_,'surprise_shocks') && ~isempty(M_.surprise_shocks)
        fprintf('dynare_estimation_init: OccBin smoother/filter: previous shocks(surprise) block detected. Deleting its content.\n')
        options_.occbin.simul.SHOCKS=zeros(1,M_.exo_nbr);
        options_.occbin.simul.exo_pos=1:M_.exo_nbr;
        M_.surprise_shocks=[];
    end
end
