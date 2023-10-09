function [fval, info, exit_flag, df, junkHessian, oo_, M_] = objective_function(xparam, BoundsInfo, oo_, estim_params_, M_, options_mom_)
% [fval, info, exit_flag, df, junk1, oo_, M_] = objective_function(xparam, BoundsInfo, oo_, estim_params_, M_, options_mom_)
% -------------------------------------------------------------------------
% This function evaluates the objective function for method of moments estimation
% =========================================================================
% INPUTS
%  o xparam:         [vector]    current value of estimated parameters as returned by set_prior()
%  o BoundsInfo:     [structure] containing parameter bounds
%  o oo_:            [structure] for results
%  o estim_params_:  [structure] describing the estimated_parameters
%  o M_              [structure] describing the model
%  o options_mom_:   [structure] information about all settings (specified by the user, preprocessor, and taken from global options_)
% -------------------------------------------------------------------------
% OUTPUTS
%  o fval:         [double] value of the quadratic form of the moment difference (except for lsqnonlin, where this is done implicitly)
%  o info:         [vector] information on error codes and penalties
%  o exit_flag:    [double] flag for exit status (0 if error, 1 if no error)
%  o df:           [matrix] analytical jacobian of the moment difference (wrt paramters), currently for GMM only
%  o junkHessian:  [matrix] empty matrix required for optimizer interface (Hessian would typically go here)
%  o oo_:          [structure] results with the following updated fields:
%                    - oo_.mom.model_moments: [vector] model moments
%                    - oo_.mom.Q: [double] value of the quadratic form of the moment difference
%                    - oo_.mom.model_moments_params_derivs: [matrix] analytical jacobian of the model moments wrt estimated parameters (currently for GMM only)
%  o M_:           [structure] updated model structure
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
% =========================================================================
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
% =========================================================================

%% TO DO
% check the info values and make use of meaningful penalties
% how do we do the penalty for the prior??


%------------------------------------------------------------------------------
% Initialization of the returned variables and others...
%------------------------------------------------------------------------------
junkHessian = [];
df = []; % required to be empty by e.g. newrat
if strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')
    if options_mom_.mom.compute_derivs && options_mom_.mom.analytic_jacobian
        if options_mom_.mom.vector_output == 1
            if options_mom_.mom.penalized_estimator
                df = nan(size(oo_.mom.data_moments,1)+length(xparam),length(xparam));
            else
                df = nan(size(oo_.mom.data_moments,1),length(xparam));
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
       fval = ones(size(oo_.mom.data_moments,1),1)*options_mom_.huge_number;
    end
    return
end


%--------------------------------------------------------------------------
% Call resol to compute steady state and model solution
%--------------------------------------------------------------------------

% Compute linear approximation around the deterministic steady state
[oo_.dr, info, M_.params] = resol(0, M_, options_mom_, oo_.dr ,oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
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
            fval = ones(size(oo_.mom.data_moments,1),1)*options_mom_.huge_number;
        end
        return
    else
        fval = Inf;
        info(4) = 0.1;
        exit_flag = 0;
        if options_mom_.mom.vector_output == 1 % lsqnonlin requires vector output
            fval = ones(size(oo_.mom.data_moments,1),1)*options_mom_.huge_number;
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
        oo_.dr.derivs = identification.get_perturbation_params_derivs(M_, options_mom_, estim_params_, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state, indpmodel, indpstderr, indpcorr, 0); %analytic derivatives of perturbation matrices
        oo_.mom.model_moments_params_derivs = NaN(options_mom_.mom.mom_nbr,totparam_nbr);
        pruned_state_space = pruned_SS.pruned_state_space_system(M_, options_mom_, oo_.dr, oo_.mom.obs_var, options_mom_.ar, 0, 1);
    else
        pruned_state_space = pruned_SS.pruned_state_space_system(M_, options_mom_, oo_.dr, oo_.mom.obs_var, options_mom_.ar, 0, 0);
    end
    oo_.mom.model_moments = NaN(options_mom_.mom.mom_nbr,1);
    for jm = 1:size(M_.matched_moments,1)
        % First moments
        if ~options_mom_.prefilter && (sum(M_.matched_moments{jm,3}) == 1)
            idx1 = (oo_.mom.obs_var == find(oo_.dr.order_var==M_.matched_moments{jm,1}) );
            oo_.mom.model_moments(jm,1) = pruned_state_space.E_y(idx1);
            if options_mom_.mom.compute_derivs && ( options_mom_.mom.analytic_standard_errors || options_mom_.mom.analytic_jacobian )
                oo_.mom.model_moments_params_derivs(jm,:) = pruned_state_space.dE_y(idx1,:);
            end
        end
        % second moments
        if (sum(M_.matched_moments{jm,3}) == 2)
            idx1 = (oo_.mom.obs_var == find(oo_.dr.order_var==M_.matched_moments{jm,1}(1)) );
            idx2 = (oo_.mom.obs_var == find(oo_.dr.order_var==M_.matched_moments{jm,1}(2)) );
            if nnz(M_.matched_moments{jm,2}) == 0
                % covariance
                if options_mom_.prefilter
                    oo_.mom.model_moments(jm,1) = pruned_state_space.Var_y(idx1,idx2);
                    if options_mom_.mom.compute_derivs && ( options_mom_.mom.analytic_standard_errors || options_mom_.mom.analytic_jacobian )
                        oo_.mom.model_moments_params_derivs(jm,:) = pruned_state_space.dVar_y(idx1,idx2,:);
                    end
                else
                    oo_.mom.model_moments(jm,1) = pruned_state_space.Var_y(idx1,idx2) + pruned_state_space.E_y(idx1)*pruned_state_space.E_y(idx2)';
                    if options_mom_.mom.compute_derivs && ( options_mom_.mom.analytic_standard_errors || options_mom_.mom.analytic_jacobian )
                        for jp=1:totparam_nbr
                            oo_.mom.model_moments_params_derivs(jm,jp) = pruned_state_space.dVar_y(idx1,idx2,jp) + pruned_state_space.dE_y(idx1,jp)*pruned_state_space.E_y(idx2)' + pruned_state_space.E_y(idx1)*pruned_state_space.dE_y(idx2,jp)';
                        end
                    end
                end
            else
                % autocovariance
                lag = -M_.matched_moments{jm,2}(2); %note that leads/lags in M_.matched_moments are transformed such that first entry is always 0 and the second is a lag
                if options_mom_.prefilter
                    oo_.mom.model_moments(jm,1) = pruned_state_space.Var_yi(idx1,idx2,lag);
                    if options_mom_.mom.compute_derivs && ( options_mom_.mom.analytic_standard_errors || options_mom_.mom.analytic_jacobian )
                        oo_.mom.model_moments_params_derivs(jm,:) = pruned_state_space.dVar_yi(idx1,idx2,lag,:);
                    end
                else
                    oo_.mom.model_moments(jm,1) = pruned_state_space.Var_yi(idx1,idx2,lag) + pruned_state_space.E_y(idx1)*pruned_state_space.E_y(idx2)';
                    if options_mom_.mom.compute_derivs && ( options_mom_.mom.analytic_standard_errors || options_mom_.mom.analytic_jacobian )
                        for jp=1:totparam_nbr
                            oo_.mom.model_moments_params_derivs(jm,jp) = vec( pruned_state_space.dVar_yi(idx1,idx2,lag,jp) + pruned_state_space.dE_y(idx1,jp)*pruned_state_space.E_y(idx2)' + pruned_state_space.E_y(idx1)*pruned_state_space.dE_y(idx2,jp)');
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
    y_sim = simult_(M_, options_mom_, oo_.dr.ys, oo_.dr, scaled_shock_series, options_mom_.order);
    % provide meaningful penalty if data is nan or inf
    if any(any(isnan(y_sim))) || any(any(isinf(y_sim)))
        if options_mom_.mom.vector_output == 1 % lsqnonlin requires vector output
            fval = Inf(size(oo_.mom.Sw,1),1);
        else
            fval = Inf;
        end
        info(1)=180;
        info(4) = 0.1;
        exit_flag = 0;
        if options_mom_.mom.vector_output == 1 % lsqnonlin requires vector output
            fval = ones(size(oo_.mom.data_moments,1),1)*options_mom_.huge_number;
        end
        return
    end
    % remove burn-in and focus on observables (note that y_sim is in declaration order)
    y_sim = y_sim(oo_.dr.order_var(oo_.mom.obs_var) , end-options_mom_.mom.long+1:end)';
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
    oo_.mom.model_moments = mom.get_data_moments(y_sim, oo_.mom.obs_var, oo_.dr.inv_order_var, M_.matched_moments, options_mom_);
end


%--------------------------------------------------------------------------
% Compute quadratic target function
%--------------------------------------------------------------------------
moments_difference = oo_.mom.data_moments - oo_.mom.model_moments;

if strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')
    residuals = sqrt(options_mom_.mom.weighting_matrix_scaling_factor)*oo_.mom.Sw*moments_difference;
    oo_.mom.Q = residuals'*residuals;
    if options_mom_.mom.vector_output == 1 % lsqnonlin requires vector output
        fval = residuals;
        if options_mom_.mom.penalized_estimator
            fval=[fval;(xparam-oo_.mom.prior.mean)./sqrt(diag(oo_.mom.prior.variance))];
        end
    else
        fval = oo_.mom.Q;
        if options_mom_.mom.penalized_estimator
            fval=fval+(xparam-oo_.mom.prior.mean)'/oo_.mom.prior.variance*(xparam-oo_.mom.prior.mean);
        end
    end
    if options_mom_.mom.compute_derivs && options_mom_.mom.analytic_jacobian
        if options_mom_.mom.penalized_estimator
            dxparam1 = eye(length(xparam));
        end
        for jp=1:length(xparam)
            dmoments_difference = - oo_.mom.model_moments_params_derivs(:,jp);
            dresiduals = sqrt(options_mom_.mom.weighting_matrix_scaling_factor)*oo_.mom.Sw*dmoments_difference;
            if options_mom_.mom.vector_output == 1 % lsqnonlin requires vector output
                if options_mom_.mom.penalized_estimator
                    df(:,jp)=[dresiduals;dxparam1(:,jp)./sqrt(diag(oo_.mom.prior.variance))];
                else
                    df(:,jp) = dresiduals;
                end
            else
                df(jp,1) = dresiduals'*residuals + residuals'*dresiduals;
                if options_mom_.mom.penalized_estimator
                    df(jp,1)=df(jp,1)+(dxparam1(:,jp))'/oo_.mom.prior.variance*(xparam-oo_.mom.prior.mean)+(xparam-oo_.mom.prior.mean)'/oo_.mom.prior.variance*(dxparam1(:,jp));
                end
            end
        end
    end
end


end % main function end

