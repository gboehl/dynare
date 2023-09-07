function DynareResults = initial_estimation_checks(objective_function,xparam1,DynareDataset,DatasetInfo,Model,EstimatedParameters,DynareOptions,BayesInfo,BoundsInfo,DynareResults)
% function DynareResults = initial_estimation_checks(objective_function,xparam1,DynareDataset,DatasetInfo,Model,EstimatedParameters,DynareOptions,BayesInfo,BoundsInfo,DynareResults)
% Checks data (complex values, ML evaluation, initial values, BK conditions,..)
%
% INPUTS
%   objective_function  [function handle] of the objective function
%   xparam1:            [vector] of parameters to be estimated
%   DynareDataset:      [dseries] object storing the dataset
%   DataSetInfo:        [structure] storing informations about the sample.
%   Model:              [structure] decribing the model
%   EstimatedParameters [structure] characterizing parameters to be estimated
%   DynareOptions       [structure] describing the options
%   BayesInfo           [structure] describing the priors
%   BoundsInfo          [structure] containing prior bounds
%   DynareResults       [structure] storing the results
%
% OUTPUTS
%    DynareResults     structure of temporary results
%
% SPECIAL REQUIREMENTS
%    none

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

%get maximum number of simultaneously observed variables for stochastic
%singularity check
maximum_number_non_missing_observations=max(sum(~isnan(DynareDataset.data(2:end,:)),2));
init_number_non_missing_observations=sum(~isnan(DynareDataset.data(1,:)),2);

if DynareOptions.heteroskedastic_filter
    if DynareOptions.order>1 
        error('initial_estimation_checks:: heteroskedastic shocks are only supported with the Kalman filter/smoother')
    end
    observations_by_period=sum(~isnan(DynareDataset.data),2);
    base_shocks = find(diag(Model.Sigma_e));
    zero_shocks = ~ismember(1:Model.exo_nbr,base_shocks);
    non_zero_shocks_by_period=repmat(length(base_shocks),size(observations_by_period));
    % check periods for which base shocks are scaled to zero
    non_zero_shocks_by_period = non_zero_shocks_by_period-sum(Model.heteroskedastic_shocks.Qscale(base_shocks,1:DynareOptions.nobs)==0,1)';
    % check periods for which base shocks are set to zero
    non_zero_shocks_by_period = non_zero_shocks_by_period-sum(Model.heteroskedastic_shocks.Qvalue(base_shocks,1:DynareOptions.nobs)==0,1)';
    % check periods for which other shocks are set to a positive number
    non_zero_shocks_by_period = non_zero_shocks_by_period+sum(Model.heteroskedastic_shocks.Qvalue(zero_shocks,1:DynareOptions.nobs)>0,1)';
end

if DynareOptions.order>1
    if any(any(isnan(DynareDataset.data)))
        error('initial_estimation_checks:: particle filtering does not support missing observations')
    end
    if DynareOptions.prefilter==1
        error('initial_estimation_checks:: particle filtering does not support the prefilter option')
    end
    if BayesInfo.with_trend
        error('initial_estimation_checks:: particle filtering does not support trends')
    end
    if DynareOptions.order>3 && DynareOptions.particle.pruning==1
        error('initial_estimation_checks:: the particle filter with order>3 does not support pruning')
    end
    if DynareOptions.particle.pruning~=DynareOptions.pruning
        warning('initial_estimation_checks:: the pruning settings differ between the particle filter and the one used for IRFs/simulations. Make sure this is intended.')
    end
end

if DynareOptions.occbin.likelihood.status || DynareOptions.occbin.smoother.status
    if DynareOptions.prefilter
        error('initial_estimation_checks:: Occbin is incompatible with the prefilter option due to the sample mean generally not corresponding to the steady state with an occasionally binding constraint.')
    end
    if ~DynareOptions.occbin.likelihood.inversion_filter && (DynareOptions.kalman_algo==2 || DynareOptions.kalman_algo==4)
        error('initial_estimation_checks:: Occbin is incompatible with the selected univariate Kalman filter.')        
    end
    if DynareOptions.fast_kalman_filter
        error('initial_estimation_checks:: Occbin is incompatible with the fast Kalman filter.')        
    end
end

if (DynareOptions.occbin.likelihood.status && DynareOptions.occbin.likelihood.inversion_filter) || (DynareOptions.occbin.smoother.status && DynareOptions.occbin.smoother.inversion_filter)
    err_index= find(diag(Model.Sigma_e)~=0);
    if length(err_index)~=length(DynareOptions.varobs)
        fprintf('initial_estimation_checks:: The IVF requires exactly as many shocks as observables.')
    end
    var_index=find(any(isnan(DynareDataset.data)));
    if ~isempty(var_index)
        fprintf('initial_estimation_checks:: The IVF requires exactly as many shocks as observables.\n')
        fprintf('initial_estimation_checks:: The data series %s contains NaN, I am therefore dropping shock %s for these time points.\n',...
            DynareOptions.varobs{var_index},Model.exo_names{DynareOptions.occbin.likelihood.IVF_shock_observable_mapping(var_index)})
    end
end

if DynareOptions.order>1 || (DynareOptions.order==1 && ~ischar(DynareOptions.mode_compute) && DynareOptions.mode_compute==11)
    if DynareOptions.order==1 && DynareOptions.mode_compute==11
        disp_string='mode_compute=11';
    else
        disp_string='particle filtering';
    end
    if Model.H==0
        error('initial_estimation_checks:: %s requires measurement error on the observables',disp_string)
    else
        if sum(diag(Model.H)>0)<length(DynareOptions.varobs)
            error('initial_estimation_checks:: %s requires as many measurement errors as observed variables',disp_string)
        else
            [~,flag]=chol(Model.H);
            if flag
                error('initial_estimation_checks:: the measurement error matrix must be positive definite')
            end
        end
    end
end

non_zero_ME=length(EstimatedParameters.H_entries_to_check_for_positive_definiteness);

print_init_check_warning=false;
if maximum_number_non_missing_observations>Model.exo_nbr+non_zero_ME
    error(['initial_estimation_checks:: Estimation can''t take place because there are less declared shocks than observed variables!'])
end
if init_number_non_missing_observations>Model.exo_nbr+non_zero_ME
    if DynareOptions.no_init_estimation_check_first_obs
        print_init_check_warning=true;
    else
        error(['initial_estimation_checks:: Estimation can''t take place because there are less declared shocks than observed variables in first period!'])
    end
end

if DynareOptions.heteroskedastic_filter
    if any(observations_by_period>(non_zero_shocks_by_period+non_zero_ME))
        error(['initial_estimation_checks:: Estimation can''t take place because too many shocks have been calibrated with a zero variance: Check heteroskedastic block and shocks calibration!'])
    end
else
    if maximum_number_non_missing_observations>length(find(diag(Model.Sigma_e)))+non_zero_ME
        error(['initial_estimation_checks:: Estimation can''t take place because too many shocks have been calibrated with a zero variance!'])
    end
end
if init_number_non_missing_observations>length(find(diag(Model.Sigma_e)))+non_zero_ME
    if DynareOptions.no_init_estimation_check_first_obs
        print_init_check_warning=true;
    else
        error(['initial_estimation_checks:: Estimation can''t take place because too many shocks have been calibrated with a zero variance in first period!'])
    end
end
if print_init_check_warning
    fprintf('ESTIMATION_CHECKS: You decided to ignore test of stochastic singularity in first_obs.\n');
    fprintf('ESTIMATION_CHECKS: If this was not done on purpose (typically when observing a stock variable [capital] in first period, on top of its flow [investment]),\n');
    fprintf('ESTIMATION_CHECKS: it may lead to a crash or provide undesired/wrong results later on!\n');
end

if (any(BayesInfo.pshape  >0 ) && DynareOptions.mh_replic) && DynareOptions.mh_nblck<1
    error(['initial_estimation_checks:: Bayesian estimation cannot be conducted with mh_nblocks=0.'])
end

% check and display warnings if steady-state solves static model (except if diffuse_filter == 1) and if steady-state changes estimated parameters
[DynareResults.steady_state] = check_steady_state_changes_parameters(Model,EstimatedParameters,DynareResults,DynareOptions, [DynareOptions.diffuse_filter==0 DynareOptions.diffuse_filter==0] );

% check and display warning if negative values of stderr or corr params are outside unit circle for Bayesian estimation
if any(BayesInfo.pshape)
    check_prior_stderr_corr(EstimatedParameters,BayesInfo);
end

% display warning if some parameters are still NaN
test_for_deep_parameters_calibration(Model);

[lnprior,~,~,info]= priordens(xparam1,BayesInfo.pshape,BayesInfo.p6,BayesInfo.p7,BayesInfo.p3,BayesInfo.p4);
if any(info)
    fprintf('The prior density evaluated at the initial values is Inf for the following parameters: %s\n',BayesInfo.name{info,1})
    error('The initial value of the prior is -Inf')
end

if isfield(Model,'filter_initial_state') && ~isempty(Model.filter_initial_state)
    state_indices=DynareResults.dr.order_var(DynareResults.dr.restrict_var_list(BayesInfo.mf0));
    for ii=1:size(state_indices,1)
        if ~isempty(Model.filter_initial_state{state_indices(ii),1})
            try
                evaluate_expression(Model.filter_initial_state{state_indices(ii),2},Model,DynareResults)
            catch
                fprintf('Unable to evaluate the expression\n %s \nfor the filter_initial_state of variable %s\n',Model.filter_initial_state{state_indices(ii),2},Model.endo_names(state_indices(ii),:))
            end
        end
    end
end

if ~isreal(DynareDataset.data)
    error('initial_estimation_checks: the data contains complex values.')
end

% Evaluate the likelihood.
ana_deriv = DynareOptions.analytic_derivation;
DynareOptions.analytic_derivation=0;
if ~isequal(DynareOptions.mode_compute,11) || ...
        (isequal(DynareOptions.mode_compute,11) && isequal(DynareOptions.order,1))
    %shut off potentially automatic switch to diffuse filter for the
    %purpose of checking stochastic singularity
    use_univariate_filters_if_singularity_is_detected_old=DynareOptions.use_univariate_filters_if_singularity_is_detected;
    DynareOptions.use_univariate_filters_if_singularity_is_detected=0;
    [fval,info] = feval(objective_function,xparam1,DynareDataset,DatasetInfo,DynareOptions,Model,EstimatedParameters,BayesInfo,BoundsInfo,DynareResults);
    if info(1)==50
        fprintf('\ninitial_estimation_checks:: The forecast error variance in the multivariate Kalman filter became singular.\n')
        fprintf('initial_estimation_checks:: This is often a sign of stochastic singularity, but can also sometimes happen by chance\n')
        fprintf('initial_estimation_checks:: for a particular combination of parameters and data realizations.\n')
        fprintf('initial_estimation_checks:: If you think the latter is the case, you should try with different initial values for the estimated parameters.\n')
        error('initial_estimation_checks:: The forecast error variance in the multivariate Kalman filter became singular.')
    end
    if info(1)==201 || info(1)==202
        message=get_error_message(info,DynareOptions);
        error('initial_estimation_checks:: %s.',message)
    end
    %reset options
    DynareOptions.use_univariate_filters_if_singularity_is_detected=use_univariate_filters_if_singularity_is_detected_old;
else
    info=0;
    fval = 0;
end
if DynareOptions.debug
    DynareResults.likelihood_at_initial_parameters=fval;
end
DynareOptions.analytic_derivation=ana_deriv;

if DynareOptions.mode_compute==13
    error('Options mode_compute=13 is only compatible with quadratic objective functions')
end


% if DynareOptions.mode_compute==5
%     if ~strcmp(func2str(objective_function),'dsge_likelihood')
%         error('Options mode_compute=5 is not compatible with non linear filters or Dsge-VAR models!')
%     end
% end
if isnan(fval)
    error('The initial value of the likelihood is NaN')
elseif imag(fval)
    error('The initial value of the likelihood is complex')
end

if info(1) > 0
    if DynareOptions.order>1
        [eigenvalues_] = check(Model,DynareOptions, DynareResults);
        if any(abs(1-abs(eigenvalues_))<abs(DynareOptions.qz_criterium-1))
            error('Your model has at least one unit root and you are using a nonlinear filter. Please set nonlinear_filter_initialization=3.')
        end
    else
        disp('Error in computing likelihood for initial parameter values')
        print_info(info, DynareOptions.noprint, DynareOptions)
    end
end

if DynareOptions.prefilter==1
    if (~DynareOptions.loglinear && any(abs(DynareResults.steady_state(BayesInfo.mfys))>1e-9)) || (DynareOptions.loglinear && any(abs(log(DynareResults.steady_state(BayesInfo.mfys)))>1e-9))
        disp(['You are trying to estimate a model with a non zero steady state for the observed endogenous'])
        disp(['variables using demeaned data!'])
        error('You should change something in your mod file...')
    end
end

if ~isequal(DynareOptions.mode_compute,11)
    disp(['Initial value of the log posterior (or likelihood): ' num2str(-fval)]);
end

if DynareOptions.mh_tune_jscale.status && (DynareOptions.mh_tune_jscale.maxiter<DynareOptions.mh_tune_jscale.stepsize)
    warning('You specified mh_tune_jscale, but the maximum number of iterations is smaller than the step size. No update will take place.')
end

if ~isempty(DynareOptions.conditional_variance_decomposition) && ~DynareOptions.moments_varendo
    disp('The conditional_variance_decomposition-option will be ignored. You need to set moments_varendo');
end

function evaluate_expression(expression,M_,oo_)
% function evaluate_expression(expression,M_,oo_)
%evaluates expressions relying on M_ and oo_ having their original names
eval(expression);
