function [fval, info, exit_flag, df, junkHessian, Q, model_moments, model_moments_params_derivs] = objective_function(xparam, data_moments, weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% [fval, info, exit_flag, df, junkHessian, Q, model_moments, model_moments_params_derivs] = objective_function(xparam, data_moments, weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
% -------------------------------------------------------------------------
% This function evaluates the objective function for a method of moments estimation
% -------------------------------------------------------------------------
% INPUTS (same ones as in dsge_likelihood.m and dsge_var_likelihood.m)
%  - xparam:               [vector]     current value of estimated parameters as returned by set_prior()
%  - data_moments:         [vector]     data with moments/irfs to match (corresponds to dataset_ in likelihood-based estimation)
%  - weighting_info:       [structure]  storing information on weighting matrices (corresponds to dataset_info in likelihood-based estimation)
%  - options_mom_:         [structure]  information about all settings (specified by the user, preprocessor, and taken from global options_)
%  - M_                    [structure]  model information
%  - estim_params_:        [structure]  information from estimated_params block
%  - bayestopt_:           [structure]  information on the prior distributions
%  - BoundsInfo:           [structure]  parameter bounds
%  - dr:                   [structure]  reduced form model
%  - endo_steady_state:    [vector]     steady state value for endogenous variables (initval)
%  - exo_steady_state:     [vector]     steady state value for exogenous variables (initval)
%  - exo_det_steady_state: [vector]     steady state value for exogenous deterministic variables (initval)
% -------------------------------------------------------------------------
% OUTPUTS
%  - fval:                         [double]  value of the quadratic form of the moment difference (except for lsqnonlin, where this is done implicitly)
%  - info:                         [vector]  information on error codes and penalties
%  - exit_flag:                    [double]  flag for exit status (0 if error, 1 if no error)
%  - df:                           [matrix]  analytical jacobian of the moment difference (wrt paramters), currently for GMM only
%  - junkHessian:                  [matrix]  empty matrix required for optimizer interface (Hessian would typically go here)
%  - Q:                            [double]  value of the quadratic form of the moment difference
%  - model_moments:                [vector]  model moments
%  - model_moments_params_derivs:  [matrix]  analytical jacobian of the model moments wrt estimated parameters (currently for GMM only)
% -------------------------------------------------------------------------
% This function is called by
%  o mom.run
%  o dynare_minimize_objective
% -------------------------------------------------------------------------
% This function calls
% o check_bounds_and_definiteness_estimation
% o get_perturbation_params_derivs
% o mom.get_data_moments
% o pruned_state_space_system
% o resol
% o set_all_parameters
% o simult_
% -------------------------------------------------------------------------

% Copyright © 2020-2023 Dynare Team
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


%------------------------------------------------------------------------------
% Initialization of the returned variables and others...
%------------------------------------------------------------------------------
model_moments_params_derivs = [];
model_moments = [];
Q = [];
junkHessian = [];
df = []; % required to be empty by e.g. newrat
if strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')
    if options_mom_.mom.compute_derivs && options_mom_.mom.analytic_jacobian
        if options_mom_.mom.vector_output == 1
            if options_mom_.mom.penalized_estimator
                df = nan(options_mom_.mom.mom_nbr+length(xparam),length(xparam));
            else
                df = nan(options_mom_.mom.mom_nbr,length(xparam));
            end
        else
            df = nan(length(xparam),1);
        end
    end
end


%--------------------------------------------------------------------------
% Get the structural parameters and define penalties
%--------------------------------------------------------------------------
% Ensure that xparam1 is a column vector; particleswarm.m requires this.
xparam = xparam(:);
M_ = set_all_parameters(xparam, estim_params_, M_);
[fval,info,exit_flag] = check_bounds_and_definiteness_estimation(xparam, M_, estim_params_, BoundsInfo);
if info(1)
    if options_mom_.mom.vector_output == 1 % lsqnonlin requires vector output
       fval = ones(options_mom_.mom.mom_nbr,1)*options_mom_.huge_number;
    end
    return
end


%--------------------------------------------------------------------------
% Call resol to compute steady state and model solution
%--------------------------------------------------------------------------
% Compute linear approximation around the deterministic steady state
[dr, info, M_.params] = resol(0, M_, options_mom_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
% Return, with endogenous penalty when possible, if resol issues an error code
if info(1)
    if info(1) == 3 || info(1) == 4 || info(1) == 5 || info(1)==6 ||info(1) == 19 ||...
            info(1) == 20 || info(1) == 21 || info(1) == 23 || info(1) == 26 || ...
            info(1) == 81 || info(1) == 84 ||  info(1) == 85 ||  info(1) == 86
        % meaningful second entry of output that can be used
        fval = Inf;
        info(4) = info(2);
        exit_flag = 0;
        if options_mom_.mom.vector_output == 1 % lsqnonlin requires vector output
            fval = ones(options_mom_.mom.mom_nbr,1)*options_mom_.huge_number;
        end
        return
    else
        fval = Inf;
        info(4) = 0.1;
        exit_flag = 0;
        if options_mom_.mom.vector_output == 1 % lsqnonlin requires vector output
            fval = ones(options_mom_.mom.mom_nbr,1)*options_mom_.huge_number;
        end
        return
    end
end


%--------------------------------------------------------------------------
% GMM: Set up pruned state-space system and compute model moments
%--------------------------------------------------------------------------
if strcmp(options_mom_.mom.mom_method,'GMM')
    if options_mom_.mom.compute_derivs && ( options_mom_.mom.analytic_standard_errors || options_mom_.mom.analytic_jacobian )
        indpmodel = []; % initialize index for model parameters
        if ~isempty(estim_params_.param_vals)
            indpmodel = estim_params_.param_vals(:,1); % values correspond to parameters declaration order, row number corresponds to order in estimated_params
        end
        indpstderr=[]; % initialize index for stderr parameters
        if ~isempty(estim_params_.var_exo)
            indpstderr = estim_params_.var_exo(:,1); % values correspond to varexo declaration order, row number corresponds to order in estimated_params
        end
        indpcorr=[]; % initialize matrix for corr paramters
        if ~isempty(estim_params_.corrx)
            indpcorr = estim_params_.corrx(:,1:2); % values correspond to varexo declaration order, row number corresponds to order in estimated_params
        end
        if estim_params_.nvn || estim_params_.ncn % nvn is number of stderr parameters and ncn is number of corr parameters of measurement innovations as declared in estimated_params
            error('Analytic computation of standard errrors does not (yet) support measurement errors.\nInstead, define them explicitly as varexo and provide measurement equations in the model definition.\nAlternatively, use numerical standard errors.')
        end
        modparam_nbr = estim_params_.np;        % number of model parameters as declared in estimated_params
        stderrparam_nbr = estim_params_.nvx;    % number of stderr parameters
        corrparam_nbr = estim_params_.ncx;      % number of corr parameters
        totparam_nbr = stderrparam_nbr+corrparam_nbr+modparam_nbr;
        dr.derivs = identification.get_perturbation_params_derivs(M_, options_mom_, estim_params_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state, indpmodel, indpstderr, indpcorr, 0); %analytic derivatives of perturbation matrices
        model_moments_params_derivs = NaN(options_mom_.mom.mom_nbr,totparam_nbr);
        pruned_state_space = pruned_SS.pruned_state_space_system(M_, options_mom_, dr, options_mom_.mom.obs_var, options_mom_.ar, 0, 1);
    else
        pruned_state_space = pruned_SS.pruned_state_space_system(M_, options_mom_, dr, options_mom_.mom.obs_var, options_mom_.ar, 0, 0);
    end
    model_moments = NaN(options_mom_.mom.mom_nbr,1);
    for jm = 1:size(M_.matched_moments,1)
        % First moments
        if ~options_mom_.prefilter && (sum(M_.matched_moments{jm,3}) == 1)
            idx1 = (options_mom_.mom.obs_var == find(dr.order_var==M_.matched_moments{jm,1}) );
            model_moments(jm,1) = pruned_state_space.E_y(idx1);
            if options_mom_.mom.compute_derivs && ( options_mom_.mom.analytic_standard_errors || options_mom_.mom.analytic_jacobian )
                model_moments_params_derivs(jm,:) = pruned_state_space.dE_y(idx1,:);
            end
        end
        % second moments
        if (sum(M_.matched_moments{jm,3}) == 2)
            idx1 = (options_mom_.mom.obs_var == find(dr.order_var==M_.matched_moments{jm,1}(1)) );
            idx2 = (options_mom_.mom.obs_var == find(dr.order_var==M_.matched_moments{jm,1}(2)) );
            if nnz(M_.matched_moments{jm,2}) == 0
                % covariance
                if options_mom_.prefilter
                    model_moments(jm,1) = pruned_state_space.Var_y(idx1,idx2);
                    if options_mom_.mom.compute_derivs && ( options_mom_.mom.analytic_standard_errors || options_mom_.mom.analytic_jacobian )
                        model_moments_params_derivs(jm,:) = pruned_state_space.dVar_y(idx1,idx2,:);
                    end
                else
                    model_moments(jm,1) = pruned_state_space.Var_y(idx1,idx2) + pruned_state_space.E_y(idx1)*pruned_state_space.E_y(idx2)';
                    if options_mom_.mom.compute_derivs && ( options_mom_.mom.analytic_standard_errors || options_mom_.mom.analytic_jacobian )
                        for jp=1:totparam_nbr
                            model_moments_params_derivs(jm,jp) = pruned_state_space.dVar_y(idx1,idx2,jp) + pruned_state_space.dE_y(idx1,jp)*pruned_state_space.E_y(idx2)' + pruned_state_space.E_y(idx1)*pruned_state_space.dE_y(idx2,jp)';
                        end
                    end
                end
            else
                % autocovariance
                lag = -M_.matched_moments{jm,2}(2); %note that leads/lags in M_.matched_moments are transformed such that first entry is always 0 and the second is a lag
                if options_mom_.prefilter
                    model_moments(jm,1) = pruned_state_space.Var_yi(idx1,idx2,lag);
                    if options_mom_.mom.compute_derivs && ( options_mom_.mom.analytic_standard_errors || options_mom_.mom.analytic_jacobian )
                        model_moments_params_derivs(jm,:) = pruned_state_space.dVar_yi(idx1,idx2,lag,:);
                    end
                else
                    model_moments(jm,1) = pruned_state_space.Var_yi(idx1,idx2,lag) + pruned_state_space.E_y(idx1)*pruned_state_space.E_y(idx2)';
                    if options_mom_.mom.compute_derivs && ( options_mom_.mom.analytic_standard_errors || options_mom_.mom.analytic_jacobian )
                        for jp=1:totparam_nbr
                            model_moments_params_derivs(jm,jp) = vec( pruned_state_space.dVar_yi(idx1,idx2,lag,jp) + pruned_state_space.dE_y(idx1,jp)*pruned_state_space.E_y(idx2)' + pruned_state_space.E_y(idx1)*pruned_state_space.dE_y(idx2,jp)');
                        end
                    end
                end
            end
        end
    end
end


%------------------------------------------------------------------------------
% SMM: Compute Moments of the model solution for Gaussian innovations
%------------------------------------------------------------------------------
if strcmp(options_mom_.mom.mom_method,'SMM')
    % create shock series with correct covariance matrix from iid standard normal shocks
    i_exo_var = setdiff(1:M_.exo_nbr, find(diag(M_.Sigma_e) == 0 )); % find singular entries in covariance
    chol_S = chol(M_.Sigma_e(i_exo_var,i_exo_var));
    scaled_shock_series = zeros(size(options_mom_.mom.shock_series)); % initialize
    scaled_shock_series(:,i_exo_var) = options_mom_.mom.shock_series(:,i_exo_var)*chol_S; % set non-zero entries
    % simulate series
    y_sim = simult_(M_, options_mom_, dr.ys, dr, scaled_shock_series, options_mom_.order);
    % provide meaningful penalty if data is nan or inf
    if any(any(isnan(y_sim))) || any(any(isinf(y_sim)))
        if options_mom_.mom.vector_output == 1 % lsqnonlin requires vector output
            fval = Inf(size(weighting_info.Sw,1),1);
        else
            fval = Inf;
        end
        info(1)=180;
        info(4) = 0.1;
        exit_flag = 0;
        if options_mom_.mom.vector_output == 1 % lsqnonlin requires vector output
            fval = ones(options_mom_.mom.mom_nbr,1)*options_mom_.huge_number;
        end
        return
    end
    % remove burn-in and focus on observables (note that y_sim is in declaration order)
    y_sim = y_sim(dr.order_var(options_mom_.mom.obs_var) , end-options_mom_.mom.long+1:end)';
    if ~all(diag(M_.H)==0)
        i_ME = setdiff(1:size(M_.H,1),find(diag(M_.H) == 0)); % find ME with 0 variance
        chol_S = chol(M_.H(i_ME,i_ME)); % decompose rest
        shock_mat=zeros(size(options_mom_.mom.ME_shock_series)); % initialize
        shock_mat(:,i_ME)=options_mom_.mom.ME_shock_series(:,i_ME)*chol_S;
        y_sim = y_sim+shock_mat;
    end
    % remove mean if centered moments
    if options_mom_.prefilter
        y_sim = bsxfun(@minus, y_sim, mean(y_sim,1));
    end
    model_moments = mom.get_data_moments(y_sim, options_mom_.mom.obs_var, dr.inv_order_var, M_.matched_moments, options_mom_);
end


%--------------------------------------------------------------------------
% Compute quadratic target function
%--------------------------------------------------------------------------
moments_difference = data_moments - model_moments;

if strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')
    residuals = sqrt(options_mom_.mom.weighting_matrix_scaling_factor)*weighting_info.Sw*moments_difference;
    Q = residuals'*residuals;
    if options_mom_.mom.vector_output == 1 % lsqnonlin requires vector output
        fval = residuals;
        if options_mom_.mom.penalized_estimator
            fval=[fval;(xparam-bayestopt_.p1)./sqrt(diag(diag(bayestopt_.p2.^2)))];
        end
    else
        fval = Q;
        if options_mom_.mom.penalized_estimator
            fval=fval+(xparam-bayestopt_.p1)'/(diag(bayestopt_.p2.^2))*(xparam-bayestopt_.p1);
        end
    end
    if options_mom_.mom.compute_derivs && options_mom_.mom.analytic_jacobian
        if options_mom_.mom.penalized_estimator
            dxparam1 = eye(length(xparam));
        end
        for jp=1:length(xparam)
            dmoments_difference = - model_moments_params_derivs(:,jp);
            dresiduals = sqrt(options_mom_.mom.weighting_matrix_scaling_factor)*weighting_info.Sw*dmoments_difference;
            if options_mom_.mom.vector_output == 1 % lsqnonlin requires vector output
                if options_mom_.mom.penalized_estimator
                    df(:,jp)=[dresiduals;dxparam1(:,jp)./sqrt(diag(diag(bayestopt_.p2.^2)))];
                else
                    df(:,jp) = dresiduals;
                end
            else
                df(jp,1) = dresiduals'*residuals + residuals'*dresiduals;
                if options_mom_.mom.penalized_estimator
                    df(jp,1)=df(jp,1)+(dxparam1(:,jp))'/(diag(bayestopt_.p2.^2))*(xparam-bayestopt_.p1)+(xparam-bayestopt_.p1)'/(diag(bayestopt_.p2.^2))*(dxparam1(:,jp));
                end
            end
        end
    end
end


end % main function end