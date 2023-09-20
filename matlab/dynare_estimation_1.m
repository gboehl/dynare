function dynare_estimation_1(var_list_,dname)
% function dynare_estimation_1(var_list_,dname)
% runs the estimation of the model
%
% INPUTS
%   var_list_:  selected endogenous variables vector
%   dname:      alternative directory name
%
% OUTPUTS
%   none
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

global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info

dispString = 'Estimation::mcmc';

if ~exist([M_.dname filesep 'Output'],'dir')
    if isoctave && octave_ver_less_than('7') && ~exist(M_.dname)
        % See https://savannah.gnu.org/bugs/index.php?61166
        % This workaround is needed for recursive estimation.
        mkdir(M_.dname)
    end
    mkdir(M_.dname,'Output');
end

if isempty(estim_params_)
    mode_compute_o = options_.mode_compute;
    mh_replic_o = options_.mh_replic;
    options_.mode_compute = 0;
    options_.mh_replic = 0;
    reset_options_related_to_estimation = true;
else
    reset_options_related_to_estimation = false;
end

%store qz_criterium
qz_criterium_old=options_.qz_criterium;
if isnan(options_.first_obs)
    first_obs_nan_indicator=true;
else
    first_obs_nan_indicator=false;
end

% Set particle filter flag.
if options_.order > 1
    if options_.particle.status
        skipline()
        disp('Estimation using a non linear filter!')
        skipline()
        if ~options_.nointeractive && ismember(options_.mode_compute,[1,3,4]) && ~strcmpi(options_.particle.filter_algorithm,'gf')% Known gradient-based optimizers
            disp('You are using a gradient-based mode-finder. Particle filtering introduces discontinuities in the')
            disp('objective function w.r.t the parameters. Thus, should use a non-gradient based optimizer.')
            fprintf('\nPlease choose a mode-finder:\n')
            fprintf('\t 0 - Continue using gradient-based method (it is most likely that you will no get any sensible result).\n')
            fprintf('\t 6 - Monte Carlo based algorithm\n')
            fprintf('\t 7 - Nelder-Mead simplex based optimization routine (Matlab optimization toolbox required)\n')
            fprintf('\t 8 - Nelder-Mead simplex based optimization routine (Dynare''s implementation)\n')
            fprintf('\t 9 - CMA-ES (Covariance Matrix Adaptation Evolution Strategy) algorithm\n')
            choice = [];
            while isempty(choice)
                choice = input('Please enter your choice: ');
                if isnumeric(choice) && isint(choice) && ismember(choice,[0 6 7 8 9])
                    if choice
                        options_.mode_compute = choice;
                    end
                else
                    fprintf('\nThis is an invalid choice (you have to choose between 0, 6, 7, 8 and 9).\n')
                    choice = [];
                end
            end
        end
    else
        error('For estimating the model with a second order approximation using a non linear filter, one should have options_.particle.status=true;')
    end
end

if ~options_.dsge_var
    if options_.particle.status
        objective_function = str2func('non_linear_dsge_likelihood');
        [options_.particle] = check_particle_filter_options(options_.particle);
    else
        if options_.occbin.likelihood.status && options_.occbin.likelihood.inversion_filter
            objective_function = str2func('occbin.IVF_posterior');
        elseif options_.conditional_likelihood.status && options_.conditional_likelihood.order==1
            objective_function = str2func('dsge_conditional_likelihood_1');
        else
            objective_function = str2func('dsge_likelihood');
        end
    end
else
    objective_function = str2func('dsge_var_likelihood');
end

[dataset_, dataset_info, xparam1, hh, M_, options_, oo_, estim_params_, bayestopt_, bounds] = ...
    dynare_estimation_init(var_list_, dname, [], M_, options_, oo_, estim_params_, bayestopt_);

if options_.dsge_var
    check_dsge_var_model(M_, estim_params_, bayestopt_);
    if dataset_info.missing.state
        error('Estimation::DsgeVarLikelihood: I cannot estimate a DSGE-VAR model with missing observations!')
    end
    if options_.noconstant
        dataset_info=var_sample_moments(options_.dsge_varlag, -1, dataset_, dataset_info);
    else
        % The steady state is non zero ==> a constant in the VAR is needed!
        dataset_info=var_sample_moments(options_.dsge_varlag, 0, dataset_, dataset_info);
    end
end

% Set sigma_e_is_diagonal flag (needed if the shocks block is not declared in the mod file).
M_.sigma_e_is_diagonal = true;
if estim_params_.ncx || any(nnz(tril(M_.Correlation_matrix,-1))) || isfield(estim_params_,'calibrated_covariances')
    M_.sigma_e_is_diagonal = false;
end

data = dataset_.data;
rawdata = dataset_info.rawdata;
data_index = dataset_info.missing.aindex;
missing_value = dataset_info.missing.state;

% Set number of observations
gend = dataset_.nobs;

% Set the number of observed variables.
n_varobs = length(options_.varobs);

% Get the number of parameters to be estimated.
nvx = estim_params_.nvx;  % Variance of the structural innovations (number of parameters).
nvn = estim_params_.nvn;  % Variance of the measurement innovations (number of parameters).
ncx = estim_params_.ncx;  % Covariance of the structural innovations (number of parameters).
ncn = estim_params_.ncn;  % Covariance of the measurement innovations (number of parameters).
np  = estim_params_.np ;  % Number of deep parameters.
nx  = nvx+nvn+ncx+ncn+np; % Total number of parameters to be estimated.

dr = oo_.dr;

if ~isempty(estim_params_)
    M_ = set_all_parameters(xparam1,estim_params_,M_);
end


%% perform initial estimation checks;
try
    oo_ = initial_estimation_checks(objective_function,xparam1,dataset_,dataset_info,M_,estim_params_,options_,bayestopt_,bounds,oo_);
catch % if check fails, provide info on using calibration if present
    e = lasterror();
    if estim_params_.full_calibration_detected %calibrated model present and no explicit starting values
        skipline(1);
        fprintf('ESTIMATION_CHECKS: There was an error in computing the likelihood for initial parameter values.\n')
        fprintf('ESTIMATION_CHECKS: If this is not a problem with the setting of options (check the error message below),\n')
        fprintf('ESTIMATION_CHECKS: you should try using the calibrated version of the model as starting values. To do\n')
        fprintf('ESTIMATION_CHECKS: this, add an empty estimated_params_init-block with use_calibration option immediately before the estimation\n')
        fprintf('ESTIMATION_CHECKS: command (and after the estimated_params-block so that it does not get overwritten):\n');
        skipline(2);
    end
    rethrow(e);
end

if isequal(options_.mode_compute,0) && isempty(options_.mode_file) && ~options_.mh_posterior_mode_estimation
    if options_.order==1 && ~options_.particle.status
        if options_.smoother
            if options_.occbin.smoother.status && options_.occbin.smoother.inversion_filter
                [~, info, ~, ~, ~, ~, ~, ~, ~, ~, oo_, atT, innov] = occbin.IVF_posterior(xparam1,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,prior_bounds(bayestopt_,options_.prior_trunc),oo_);
                if ismember(info(1),[303,304,306])
                    fprintf('\nIVF: smoother did not succeed. No results will be written to oo_.\n')
                else
                    updated_variables = atT*nan;
                    measurement_error=[];
                    ys = oo_.dr.ys;
                    trend_coeff = zeros(length(options_.varobs_id),1);
                    bayestopt_.mf = bayestopt_.smoother_var_list(bayestopt_.smoother_mf);
                    options_nk=options_.nk;
                    options_.nk=[]; %unset options_.nk and reset it later
                    [oo_]=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,atT,innov,measurement_error,updated_variables,ys,trend_coeff);
                    options_.nk=options_nk;
                end
            else
                if options_.occbin.smoother.status
                    [atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,T,R,P,PK,decomp,Trend,state_uncertainty,M_,oo_,bayestopt_] = occbin.DSGE_smoother(xparam1,gend,transpose(data),data_index,missing_value,M_,oo_,options_,bayestopt_,estim_params_,dataset_,dataset_info);
                    if oo_.occbin.smoother.error_flag(1)==0
                        [oo_]=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,P,PK,decomp,Trend,state_uncertainty);
                    else
                        fprintf('\nOccbin: smoother did not succeed. No results will be written to oo_.\n')
                    end
                else
                    [atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,T,R,P,PK,decomp,Trend,state_uncertainty,M_,oo_,bayestopt_] = DsgeSmoother(xparam1,gend,transpose(data),data_index,missing_value,M_,oo_,options_,bayestopt_,estim_params_);
                    [oo_]=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,P,PK,decomp,Trend,state_uncertainty);
                end
            end
            if options_.forecast > 0
                oo_.forecast = dyn_forecast(var_list_,M_,options_,oo_,'smoother',dataset_info);
            end
        end
        %reset qz_criterium
        options_.qz_criterium=qz_criterium_old;
        return
    else %allow to continue, e.g. with MCMC_jumping_covariance
        if options_.smoother
            error('Estimation:: Particle Smoothers are not yet implemented.')
        end
    end
end

%% Estimation of the posterior mode or likelihood mode



if ~isequal(options_.mode_compute,0) && ~options_.mh_posterior_mode_estimation
    optimizer_vec = [options_.mode_compute;num2cell(options_.additional_optimizer_steps)];
    for optim_iter = 1:length(optimizer_vec)
        current_optimizer = optimizer_vec{optim_iter};

        [xparam1, fval, exitflag, hh, options_, Scale, new_rat_hess_info] = dynare_minimize_objective(objective_function,xparam1,current_optimizer,options_,[bounds.lb bounds.ub],bayestopt_.name,bayestopt_,hh,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_);
        fprintf('\nFinal value of minus the log posterior (or likelihood):%f \n', fval);

        if isnumeric(current_optimizer)
            if current_optimizer==5
                newratflag = new_rat_hess_info.newratflag;
                new_rat_hess_info = new_rat_hess_info.new_rat_hess_info;
            elseif current_optimizer==6 %save scaling factor
                save([M_.dname filesep 'Output' filesep M_.fname '_optimal_mh_scale_parameter.mat'],'Scale');
                options_.mh_jscale = Scale;
                bayestopt_.jscale(:) = options_.mh_jscale;
            end
        end
        if ~isnumeric(current_optimizer) || ~isequal(current_optimizer,6) %always already computes covariance matrix
            if options_.cova_compute == 1 %user did not request covariance not to be computed
                if options_.analytic_derivation && strcmp(func2str(objective_function),'dsge_likelihood')
                    ana_deriv_old = options_.analytic_derivation;
                    options_.analytic_derivation = 2;
                    [~,~,~,~,hh] = feval(objective_function,xparam1, ...
                        dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_);
                    options_.analytic_derivation = ana_deriv_old;
                elseif ~isnumeric(current_optimizer) || ~(isequal(current_optimizer,5) && newratflag~=1 && strcmp(func2str(objective_function),'dsge_likelihood'))
                    % enter here if i) not mode_compute_5, ii) if mode_compute_5 and newratflag==1;
                    % with flag==0 or 2 and dsge_likelihood, we force to use
                    % the hessian from outer product gradient of optimizer 5 below
                    if options_.hessian.use_penalized_objective
                        penalized_objective_function = str2func('penalty_objective_function');
                        hh = hessian(penalized_objective_function, xparam1, options_.gstep, objective_function, fval, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, bounds,oo_);
                    else
                        hh = hessian(objective_function, xparam1, options_.gstep, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, bounds,oo_);
                    end
                    hh = reshape(hh, nx, nx);
                elseif isnumeric(current_optimizer) && isequal(current_optimizer,5)
                    % other numerical hessian options available with optimizer
                    % 5 and dsge_likelihood
                    %
                    % if newratflag == 0
                    % compute outer product gradient of optimizer 5
                    %
                    % if newratflag == 2
                    % compute 'mixed' outer product gradient of optimizer 5
                    % with diagonal elements computed with numerical second order derivatives
                    %
                    % uses univariate filters, so to get max # of available
                    % densities for outer product gradient
                    kalman_algo0 = options_.kalman_algo;
                    compute_hessian = 1;
                    if ~((options_.kalman_algo == 2) || (options_.kalman_algo == 4))
                        options_.kalman_algo=2;
                        if options_.lik_init == 3
                            options_.kalman_algo=4;
                        end
                    elseif newratflag==0 % hh already contains outer product gradient with univariate filter
                        compute_hessian = 0;
                    end
                    if compute_hessian
                        crit = options_.newrat.tolerance.f;
                        newratflag = newratflag>0;
                        hh = reshape(mr_hessian(xparam1,objective_function,fval,newratflag,crit,new_rat_hess_info,[bounds.lb bounds.ub],bayestopt_.p2,0,dataset_, dataset_info, options_,M_,estim_params_,bayestopt_,bounds,oo_), nx, nx);
                    end
                    options_.kalman_algo = kalman_algo0;
                end
            end
        end
        parameter_names = bayestopt_.name;
    end
    if options_.cova_compute || current_optimizer==5 || current_optimizer==6
        save([M_.dname filesep 'Output' filesep M_.fname '_mode.mat'],'xparam1','hh','parameter_names','fval');
    else
        save([M_.dname filesep 'Output' filesep M_.fname '_mode.mat'],'xparam1','parameter_names','fval');
    end
end

if ~options_.mh_posterior_mode_estimation && options_.cova_compute
    check_hessian_at_the_mode(hh, xparam1, M_, estim_params_, options_, bounds);
end

if options_.mode_check.status && ~options_.mh_posterior_mode_estimation
    ana_deriv_old = options_.analytic_derivation;
    options_.analytic_derivation = 0;
    mode_check(objective_function,xparam1,hh,options_,M_,estim_params_,bayestopt_,bounds,false,...
               dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, bounds,oo_);
    options_.analytic_derivation = ana_deriv_old;
end

oo_.posterior.optimization.mode = [];
oo_.posterior.optimization.Variance = [];
oo_.posterior.optimization.log_density=[];

invhess=[];
if ~options_.mh_posterior_mode_estimation
    oo_.posterior.optimization.mode = xparam1;
    if exist('fval','var')
        oo_.posterior.optimization.log_density=-fval;
    end
    if options_.cova_compute
        hsd = sqrt(diag(hh));
        invhess = inv(hh./(hsd*hsd'))./(hsd*hsd');
        stdh = sqrt(diag(invhess));
        oo_.posterior.optimization.Variance = invhess;
    end
else
    variances = bayestopt_.p2.*bayestopt_.p2;
    idInf = isinf(variances);
    variances(idInf) = 1;
    invhess = options_.mh_posterior_mode_estimation*diag(variances);
    xparam1 = bayestopt_.p5;
    idNaN = isnan(xparam1);
    xparam1(idNaN) = bayestopt_.p1(idNaN);
    outside_bound_pars=find(xparam1 < bounds.lb | xparam1 > bounds.ub);
    xparam1(outside_bound_pars) = bayestopt_.p1(outside_bound_pars);
end

if ~options_.cova_compute
    stdh = NaN(length(xparam1),1);
end

if options_.particle.status && isfield(options_.particle,'posterior_sampler')
    if strcmpi(options_.particle.posterior_sampler,'Herbst_Schorfheide')
        Herbst_Schorfheide_sampler(objective_function,xparam1,bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_)
    elseif strcmpi(options_.particle.posterior_sampler,'DSMH')
        DSMH_sampler(objective_function,xparam1,bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_)
    end
end
        
if any(bayestopt_.pshape > 0) && ~options_.mh_posterior_mode_estimation
    % display results table and store parameter estimates and standard errors in results
    oo_ = display_estimation_results_table(xparam1, stdh, M_, options_, estim_params_, bayestopt_, oo_, prior_dist_names, 'Posterior', 'posterior');
    % Laplace approximation to the marginal log density:
    if options_.cova_compute
        estim_params_nbr = size(xparam1,1);
        if ispd(invhess)
            log_det_invhess = log(det(invhess./(stdh*stdh')))+2*sum(log(stdh));
            likelihood = feval(objective_function,xparam1,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_);
            oo_.MarginalDensity.LaplaceApproximation = .5*estim_params_nbr*log(2*pi) + .5*log_det_invhess - likelihood;
        else
            oo_.MarginalDensity.LaplaceApproximation = NaN;
        end
        skipline()
        disp(sprintf('Log data density [Laplace approximation] is %f.',oo_.MarginalDensity.LaplaceApproximation))
        skipline()
    end
    if options_.dsge_var
        [~,~,~,~,~,~,~,oo_.dsge_var.posterior_mode.PHI_tilde,oo_.dsge_var.posterior_mode.SIGMA_u_tilde,oo_.dsge_var.posterior_mode.iXX,oo_.dsge_var.posterior_mode.prior] =...
            feval(objective_function,xparam1,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_);
    end

elseif ~any(bayestopt_.pshape > 0) && ~options_.mh_posterior_mode_estimation
    oo_=display_estimation_results_table(xparam1, stdh, M_, options_, estim_params_, bayestopt_, oo_, prior_dist_names, 'Maximum Likelihood', 'mle');
end

if np > 0
    pindx = estim_params_.param_vals(:,1);
    save([M_.dname filesep 'Output' filesep M_.fname '_params.mat'],'pindx');
end

invhess = set_mcmc_jumping_covariance(invhess, nx, options_.MCMC_jumping_covariance, bayestopt_, 'dynare_estimation_1');

if (any(bayestopt_.pshape  >0 ) && options_.mh_replic) || ...
        (any(bayestopt_.pshape >0 ) && options_.load_mh_file)  %% not ML estimation
    %reset bounds as lb and ub must only be operational during mode-finding
    bounds = set_mcmc_prior_bounds(xparam1, bayestopt_, options_, 'dynare_estimation_1');
    % Tunes the jumping distribution's scale parameter
    if options_.mh_tune_jscale.status
        if strcmp(options_.posterior_sampler_options.posterior_sampling_method, 'random_walk_metropolis_hastings')
            options_.mh_jscale = tune_mcmc_mh_jscale_wrapper(invhess, options_, M_, objective_function, xparam1, bounds,...
                                                             dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, bounds, oo_);
            bayestopt_.jscale(:) = options_.mh_jscale;
            fprintf('mh_jscale has been set equal to %s\n', num2str(options_.mh_jscale));
        else
            warning('mh_tune_jscale is only available with Random Walk Metropolis Hastings!')
        end
    end
    % runs MCMC
    if options_.mh_replic || options_.load_mh_file
        posterior_sampler_options = options_.posterior_sampler_options.current_options;
        posterior_sampler_options.invhess = invhess;
        [posterior_sampler_options, options_, bayestopt_] = check_posterior_sampler_options(posterior_sampler_options, M_.fname, M_.dname, options_, bounds, bayestopt_);
        % store current options in global
        options_.posterior_sampler_options.current_options = posterior_sampler_options;
        if options_.mh_replic
            ana_deriv_old = options_.analytic_derivation;
            options_.analytic_derivation = 0;
            posterior_sampler(objective_function,posterior_sampler_options.proposal_distribution,xparam1,posterior_sampler_options,bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_,dispString);
            options_.analytic_derivation = ana_deriv_old;
        end
    end
    %% Here I discard first mh_drop percent of the draws:
    CutSample(M_, options_, dispString);
    if options_.mh_posterior_mode_estimation
        [~,~,posterior_mode,~] = compute_mh_covariance_matrix(bayestopt_,M_.fname,M_.dname);
        oo_=fill_mh_mode(posterior_mode',NaN(length(posterior_mode),1),M_,options_,estim_params_,bayestopt_,oo_,'posterior');
        %reset qz_criterium
        options_.qz_criterium=qz_criterium_old;
        return
    else
        %get stored results if required
        if options_.load_mh_file && options_.load_results_after_load_mh
            oo_load_mh=load([M_.dname filesep 'Output' filesep M_.fname '_results'],'oo_');
        end
        if ~options_.nodiagnostic
            if (options_.mh_replic>0 || (options_.load_mh_file && ~options_.load_results_after_load_mh))
                oo_= mcmc_diagnostics(options_, estim_params_, M_,oo_);
            elseif options_.load_mh_file && options_.load_results_after_load_mh
                if isfield(oo_load_mh.oo_,'convergence')
                    oo_.convergence=oo_load_mh.oo_.convergence;
                end
            end
        end
        %% Estimation of the marginal density from the Mh draws:
        if options_.mh_replic || (options_.load_mh_file && ~options_.load_results_after_load_mh)
            [marginal,oo_] = marginal_density(M_, options_, estim_params_, oo_, bayestopt_);
            % Store posterior statistics by parameter name
            oo_ = GetPosteriorParametersStatistics(estim_params_, M_, options_, bayestopt_, oo_, prior_dist_names);
            if ~options_.nograph
                oo_ = PlotPosteriorDistributions(estim_params_, M_, options_, bayestopt_, oo_);
            end
            % Store posterior mean in a vector and posterior variance in
            % a matrix
            [oo_.posterior.metropolis.mean,oo_.posterior.metropolis.Variance] ...
                = GetPosteriorMeanVariance(M_,options_.mh_drop);
        elseif options_.load_mh_file && options_.load_results_after_load_mh
            %% load fields from previous MCMC run stored in results-file
            field_names={'posterior_mode','posterior_std_at_mode',...% fields set by marginal_density
                         'posterior_mean','posterior_hpdinf','posterior_hpdsup','posterior_median','posterior_variance','posterior_std','posterior_deciles','posterior_density',...% fields set by GetPosteriorParametersStatistics
                         'prior_density',...%fields set by PlotPosteriorDistributions
                        };
            for field_iter=1:size(field_names,2)
                if isfield(oo_load_mh.oo_,field_names{1,field_iter})
                    oo_.(field_names{1,field_iter})=oo_load_mh.oo_.(field_names{1,field_iter});
                end
            end
            % field set by marginal_density
            if isfield(oo_load_mh.oo_,'MarginalDensity') && isfield(oo_load_mh.oo_.MarginalDensity,'ModifiedHarmonicMean')
                oo_.MarginalDensity.ModifiedHarmonicMean=oo_load_mh.oo_.MarginalDensity.ModifiedHarmonicMean;
            end
            % field set by GetPosteriorMeanVariance
            if isfield(oo_load_mh.oo_,'posterior') && isfield(oo_load_mh.oo_.posterior,'metropolis')
                oo_.posterior.metropolis=oo_load_mh.oo_.posterior.metropolis;
            end
        end
        [error_flag,~,options_]= metropolis_draw(1,options_,estim_params_,M_);
        if ~(~isempty(options_.sub_draws) && options_.sub_draws==0)
            if options_.bayesian_irf
                if error_flag
                    error('%s: I cannot compute the posterior IRFs!',dispString)
                end
                PosteriorIRF('posterior',dispString);
            end
            if options_.moments_varendo
                if error_flag
                    error('%s: I cannot compute the posterior moments for the endogenous variables!',dispString)
                end
                if options_.load_mh_file && options_.mh_replic==0 %user wants to recompute results
                   [MetropolisFolder, info] = CheckPath('metropolis',M_.dname);
                   if ~info
                       generic_post_data_file_name={'Posterior2ndOrderMoments','decomposition','PosteriorVarianceDecomposition','correlation','PosteriorCorrelations','conditional decomposition','PosteriorConditionalVarianceDecomposition'};
                       for ii=1:length(generic_post_data_file_name)
                           delete_stale_file([MetropolisFolder filesep M_.fname '_' generic_post_data_file_name{1,ii} '*']);
                       end
                       % restore compatibility for loading pre-4.6.2 runs where estim_params_ was not saved; see 6e06acc7 and !1944
                       NumberOfDrawsFiles = length(dir([M_.dname '/metropolis/' M_.fname '_posterior_draws*' ]));
                       if NumberOfDrawsFiles>0
                           temp=load([M_.dname '/metropolis/' M_.fname '_posterior_draws1']);
                           if ~isfield(temp,'estim_params_')
                               for file_iter=1:NumberOfDrawsFiles
                                   save([M_.dname '/metropolis/' M_.fname '_posterior_draws' num2str(file_iter)],'estim_params_','-append') 
                               end
                           end
                       end
                   end
                end
                oo_ = compute_moments_varendo('posterior',options_,M_,oo_,var_list_);
            end
            if options_.smoother || ~isempty(options_.filter_step_ahead) || options_.forecast
                if error_flag
                    error('%s: I cannot compute the posterior statistics!',dispString)
                end
                if options_.order==1 && ~options_.particle.status
                    prior_posterior_statistics('posterior',dataset_,dataset_info,dispString); %get smoothed and filtered objects and forecasts
                else
                    error('%s: Particle Smoothers are not yet implemented.',dispString)
                end
            end
        else
            fprintf('%s: sub_draws was set to 0. Skipping posterior computations.',dispString);
        end
        xparam1 = get_posterior_parameters('mean',M_,estim_params_,oo_,options_);
        M_ = set_all_parameters(xparam1,estim_params_,M_);
    end
end

if options_.particle.status
    %reset qz_criterium
    options_.qz_criterium=qz_criterium_old;
    return
end

if (~((any(bayestopt_.pshape > 0) && options_.mh_replic) || (any(bayestopt_.pshape> 0) && options_.load_mh_file)) ...
    || ~options_.smoother ) && ~options_.partial_information  % to be fixed
    %% ML estimation, or posterior mode without Metropolis-Hastings or Metropolis without Bayesian smoothed variables
    if options_.occbin.smoother.status && options_.occbin.smoother.inversion_filter
        [~, info, ~, ~, ~, ~, ~, ~, ~, ~, oo_, atT, innov] = occbin.IVF_posterior(xparam1,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,prior_bounds(bayestopt_,options_.prior_trunc),oo_);
        if ismember(info(1),[303,304,306])
            fprintf('\nIVF: smoother did not succeed. No results will be written to oo_.\n')
        else
            updated_variables = atT*nan;
            measurement_error=[];
            ys = oo_.dr.ys;
            trend_coeff = zeros(length(options_.varobs_id),1);
            bayestopt_.mf = bayestopt_.smoother_var_list(bayestopt_.smoother_mf);
            options_nk=options_.nk;
            options_.nk=[]; %unset options_.nk and reset it later
            [oo_, yf]=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,atT,innov,measurement_error,updated_variables,ys,trend_coeff);
            options_.nk=options_nk;
        end        
    else
        if options_.occbin.smoother.status
            [atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,T,R,P,PK,decomp,Trend,state_uncertainty,M_,oo_,bayestopt_] = occbin.DSGE_smoother(xparam1,dataset_.nobs,transpose(dataset_.data),dataset_info.missing.aindex,dataset_info.missing.state,M_,oo_,options_,bayestopt_,estim_params_,dataset_,dataset_info);
            if oo_.occbin.smoother.error_flag(1)==0
                [oo_,yf]=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,P,PK,decomp,Trend,state_uncertainty);
            else 
                fprintf('\nOccbin: smoother did not succeed. No results will be written to oo_.\n')
            end
        else
            [atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,T,R,P,PK,decomp,Trend,state_uncertainty,M_,oo_,bayestopt_] = DsgeSmoother(xparam1,dataset_.nobs,transpose(dataset_.data),dataset_info.missing.aindex,dataset_info.missing.state,M_,oo_,options_,bayestopt_,estim_params_);
            [oo_,yf]=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,P,PK,decomp,Trend,state_uncertainty);
        end
        [oo_,yf]=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,P,PK,decomp,Trend,state_uncertainty);
    end
    if ~options_.nograph
        [nbplt,nr,nc,lr,lc,nstar] = pltorg(M_.exo_nbr);
        if ~exist([M_.dname '/graphs'],'dir')
            mkdir(M_.dname,'graphs');
        end
        if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
            fidTeX = fopen([M_.dname, '/graphs/' M_.fname '_SmoothedShocks.tex'],'w');
            fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation_1.m (Dynare).\n');
            fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
            fprintf(fidTeX,' \n');
        end
        for plt = 1:nbplt
            fh = dyn_figure(options_.nodisplay,'Name','Smoothed shocks');
            NAMES = [];
            if options_.TeX, TeXNAMES = []; end
            nstar0=min(nstar,M_.exo_nbr-(plt-1)*nstar);
            if gend==1
                marker_string{1,1}='-ro';
                marker_string{2,1}='-ko';
            else
                marker_string{1,1}='-r';
                marker_string{2,1}='-k';
            end
            for i=1:nstar0
                k = (plt-1)*nstar+i;
                subplot(nr,nc,i);
                plot([1 gend],[0 0],marker_string{1,1},'linewidth',.5)
                hold on
                plot(1:gend,innov(k,:),marker_string{2,1},'linewidth',1)
                hold off
                name = M_.exo_names{k};
                if ~isempty(options_.XTick)
                    set(gca,'XTick',options_.XTick)
                    set(gca,'XTickLabel',options_.XTickLabel)
                end
                if gend>1
                    xlim([1 gend])
                end
                if options_.TeX
                    title(['$' M_.exo_names_tex{k} '$'],'Interpreter','latex')
                else
                    title(name,'Interpreter','none')
                end
            end
            dyn_saveas(fh,[M_.dname, '/graphs/' M_.fname '_SmoothedShocks' int2str(plt)],options_.nodisplay,options_.graph_format);
            if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s_SmoothedShocks%s}\n',options_.figures.textwidth*min(i/nc,1),[M_.dname, '/graphs/' M_.fname],int2str(plt));
                fprintf(fidTeX,'\\caption{Smoothed shocks.}');
                fprintf(fidTeX,'\\label{Fig:SmoothedShocks:%s}\n',int2str(plt));
                fprintf(fidTeX,'\\end{figure}\n');
                fprintf(fidTeX,'\n');
            end
        end
        if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
            fprintf(fidTeX,'\n');
            fprintf(fidTeX,'%% End of TeX file.\n');
            fclose(fidTeX);
        end
    end
    if nvn
        number_of_plots_to_draw = 0;
        index = [];
        for obs_iter=1:n_varobs
            if max(abs(measurement_error(obs_iter,:))) > options_.ME_plot_tol;
                number_of_plots_to_draw = number_of_plots_to_draw + 1;
                index = cat(1,index,obs_iter);
            end
        end
        if ~options_.nograph
            [nbplt,nr,nc,lr,lc,nstar] = pltorg(number_of_plots_to_draw);
            if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                fidTeX = fopen([M_.dname, '/graphs/' M_.fname '_SmoothedObservationErrors.tex'],'w');
                fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation_1.m (Dynare).\n');
                fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
                fprintf(fidTeX,' \n');
            end
            for plt = 1:nbplt
                fh = dyn_figure(options_.nodisplay,'Name','Smoothed observation errors');
                nstar0=min(nstar,number_of_plots_to_draw-(plt-1)*nstar);
                if gend==1
                    marker_string{1,1}='-ro';
                    marker_string{2,1}='-ko';
                else
                    marker_string{1,1}='-r';
                    marker_string{2,1}='-k';
                end
                for i=1:nstar0
                    k = (plt-1)*nstar+i;
                    subplot(nr,nc,i);
                    plot([1 gend],[0 0],marker_string{1,1},'linewidth',.5)
                    hold on
                    plot(1:gend,measurement_error(index(k),:),marker_string{2,1},'linewidth',1)
                    hold off
                    name = options_.varobs{index(k)};
                    if gend>1
                        xlim([1 gend])
                    end
                    if ~isempty(options_.XTick)
                        set(gca,'XTick',options_.XTick)
                        set(gca,'XTickLabel',options_.XTickLabel)
                    end
                    if options_.TeX
                        idx = strmatch(options_.varobs{index(k)}, M_.endo_names, 'exact');
                        title(['$' M_.endo_names_tex{idx} '$'],'Interpreter','latex')
                    else
                        title(name,'Interpreter','none')
                    end
                end
                dyn_saveas(fh,[M_.dname, '/graphs/' M_.fname '_SmoothedObservationErrors' int2str(plt)],options_.nodisplay,options_.graph_format);
                if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                    fprintf(fidTeX,'\\begin{figure}[H]\n');
                    fprintf(fidTeX,'\\centering \n');
                    fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s_SmoothedObservationErrors%s}\n',options_.figures.textwidth*min(i/nc,1),[M_.dname, '/graphs/' M_.fname],int2str(plt));
                    fprintf(fidTeX,'\\caption{Smoothed observation errors.}');
                    fprintf(fidTeX,'\\label{Fig:SmoothedObservationErrors:%s}\n',int2str(plt));
                    fprintf(fidTeX,'\\end{figure}\n');
                    fprintf(fidTeX,'\n');
                end
            end
            if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                fprintf(fidTeX,'\n');
                fprintf(fidTeX,'%% End of TeX file.\n');
                fclose(fidTeX);
            end
        end
    end
    %%
    %%  Historical and smoothed variabes
    %%
    if ~options_.nograph
        [nbplt,nr,nc,lr,lc,nstar] = pltorg(n_varobs);
        if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
            fidTeX = fopen([M_.dname, '/graphs/' M_.fname '_HistoricalAndSmoothedVariables.tex'],'w');
            fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation_1.m (Dynare).\n');
            fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
            fprintf(fidTeX,' \n');
        end
        for plt = 1:nbplt
            fh = dyn_figure(options_.nodisplay,'Name','Historical and smoothed variables');
            NAMES = [];
            nstar0=min(nstar,n_varobs-(plt-1)*nstar);
            if gend==1
                marker_string{1,1}='-ro';
                marker_string{2,1}='--ko';
            else
                marker_string{1,1}='-r';
                marker_string{2,1}='--k';
            end
            for i=1:nstar0
                k = (plt-1)*nstar+i;
                subplot(nr,nc,i);
                plot(1:gend,yf(k,:),marker_string{1,1},'linewidth',1)
                hold on
                plot(1:gend,rawdata(:,k),marker_string{2,1},'linewidth',1)
                hold off
                name = options_.varobs{k};
                if ~isempty(options_.XTick)
                    set(gca,'XTick',options_.XTick)
                    set(gca,'XTickLabel',options_.XTickLabel)
                end
                if gend>1
                    xlim([1 gend])
                end
                if options_.TeX
                    idx = strmatch(options_.varobs{k}, M_.endo_names,'exact');
                    title(['$' M_.endo_names_tex{idx} '$'],'Interpreter','latex')
                else
                    title(name,'Interpreter','none')
                end
            end
            dyn_saveas(fh,[M_.dname, '/graphs/' M_.fname '_HistoricalAndSmoothedVariables' int2str(plt)],options_.nodisplay,options_.graph_format);
            if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s_HistoricalAndSmoothedVariables%s}\n',options_.figures.textwidth*min(i/nc,1),[M_.dname, '/graphs/' M_.fname],int2str(plt));
                fprintf(fidTeX,'\\caption{Historical and smoothed variables.}');
                fprintf(fidTeX,'\\label{Fig:HistoricalAndSmoothedVariables:%s}\n',int2str(plt));
                fprintf(fidTeX,'\\end{figure}\n');
                fprintf(fidTeX,'\n');
            end
        end
        if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
            fprintf(fidTeX,'\n');
            fprintf(fidTeX,'%% End of TeX file.\n');
            fclose(fidTeX);
        end
    end
end

if options_.forecast > 0 && options_.mh_replic == 0 && ~options_.load_mh_file
    oo_.forecast = dyn_forecast(var_list_,M_,options_,oo_,'smoother',dataset_info);
end

if np > 0
    pindx = estim_params_.param_vals(:,1);
    save([M_.dname filesep 'Output' filesep M_.fname '_pindx.mat'] ,'pindx');
end

%reset qz_criterium
options_.qz_criterium=qz_criterium_old;

if reset_options_related_to_estimation
    options_.mode_compute = mode_compute_o;
    options_.mh_replic = mh_replic_o;
end
if first_obs_nan_indicator
    options_.first_obs=NaN;
end
