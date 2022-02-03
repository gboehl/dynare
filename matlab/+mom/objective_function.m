function [fval, info, exit_flag, df, junk1, oo_, M_, options_mom_] = objective_function(xparam1, Bounds, oo_, estim_params_, M_, options_mom_)
% [fval, info, exit_flag, df, junk1, oo_, M_, options_mom_] = objective_function(xparam1, Bounds, oo_, estim_params_, M_, options_mom_)
% -------------------------------------------------------------------------
% This function evaluates the objective function for GMM/SMM estimation
% =========================================================================
% INPUTS
%   o xparam1:                  current value of estimated parameters as returned by set_prior()
%   o Bounds:                   structure containing parameter bounds
%   o oo_:                      structure for results
%   o estim_params_:            structure describing the estimated_parameters
%   o M_                        structure describing the model
%   o options_mom_:             structure information about all settings (specified by the user, preprocessor, and taken from global options_)
% -------------------------------------------------------------------------
% OUTPUTS
%   o fval:                     value of the quadratic form of the moment difference (except for lsqnonlin, where this is done implicitly)
%   o info:                     vector storing error code and penalty 
%   o exit_flag:                0 if error, 1 if no error
%   o df:                       analytical parameter Jacobian of the quadratic form of the moment difference (for GMM only)
%   o junk1:                    empty matrix required for optimizer interface (Hessian would go here)
%   o oo_:                      structure containing the results with the following updated fields:
%      - mom.model_moments       [numMom x 1] vector with model moments
%      - mom.Q                   value of the quadratic form of the moment difference
%      - mom.model_moments_params_derivs
%                                [numMom x numParams] Jacobian matrix of derivatives of model_moments with respect to estimated parameters
%                                (only for GMM with analytical derivatives)
%   o M_:                       Matlab's structure describing the model
%   o options_mom_:             structure information about all settings (specified by the user, preprocessor, and taken from global options_)
% -------------------------------------------------------------------------
% This function is called by
%  o mom.run.m
%  o dynare_minimize_objective.m
% -------------------------------------------------------------------------
% This function calls
%  o check_bounds_and_definiteness_estimation
%  o pruned_state_space_system
%  o resol
%  o set_all_parameters
% =========================================================================
% Copyright (C) 2020-2021 Dynare Team
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
% -------------------------------------------------------------------------
% Author(s): 
% o Willi Mutschler (willi@mutschler.eu)
% o Johannes Pfeifer (jpfeifer@uni-koeln.de)
% =========================================================================

%------------------------------------------------------------------------------
% 0. Initialization of the returned variables and others...
%------------------------------------------------------------------------------
if options_mom_.mom.compute_derivs && options_mom_.mom.analytic_jacobian
    if options_mom_.vector_output == 1
        if options_mom_.mom.penalized_estimator
            df           = nan(size(oo_.mom.data_moments,1)+length(xparam1),length(xparam1));
        else
            df           = nan(size(oo_.mom.data_moments,1),length(xparam1));
        end
    else
        df           = nan(1,length(xparam1));
    end
else
    df=[]; %required to be empty by e.g. newrat
end
junk1        = [];
junk2        = [];

%--------------------------------------------------------------------------
% 1. Get the structural parameters & define penalties
%--------------------------------------------------------------------------

% Ensure that xparam1 is a column vector; particleswarm.m requires this.
xparam1 = xparam1(:);

M_ = set_all_parameters(xparam1, estim_params_, M_);

[fval,info,exit_flag]=check_bounds_and_definiteness_estimation(xparam1, M_, estim_params_, Bounds);
if info(1)
    if options_mom_.vector_output == 1 % lsqnonlin requires vector output
       fval = ones(size(oo_.mom.data_moments,1),1)*options_mom_.huge_number;
    end
    return
end

%--------------------------------------------------------------------------
% 2. call resol to compute steady state and model solution
%--------------------------------------------------------------------------

% Compute linear approximation around the deterministic steady state
[dr, info, M_, oo_] = resol(0, M_, options_mom_, oo_);

% Return, with endogenous penalty when possible, if resol issues an error code
if info(1)
    if info(1) == 3 || info(1) == 4 || info(1) == 5 || info(1)==6 ||info(1) == 19 ||...
            info(1) == 20 || info(1) == 21 || info(1) == 23 || info(1) == 26 || ...
            info(1) == 81 || info(1) == 84 ||  info(1) == 85 ||  info(1) == 86
        %meaningful second entry of output that can be used
        fval = Inf;
        info(4) = info(2);
        exit_flag = 0;
        if options_mom_.vector_output == 1 % lsqnonlin requires vector output
            fval = ones(size(oo_.mom.data_moments,1),1)*options_mom_.huge_number;
        end
        return
    else
        fval = Inf;
        info(4) = 0.1;
        exit_flag = 0;
        if options_mom_.vector_output == 1 % lsqnonlin requires vector output
            fval = ones(size(oo_.mom.data_moments,1),1)*options_mom_.huge_number;
        end
        return
    end
end

if strcmp(options_mom_.mom.mom_method,'GMM')
    %--------------------------------------------------------------------------
    % 3. Set up pruned state-space system and compute model moments
    %--------------------------------------------------------------------------
    if options_mom_.mom.compute_derivs && ( options_mom_.mom.analytic_standard_errors || options_mom_.mom.analytic_jacobian )
        indpmodel = []; %initialize index for model parameters
        if ~isempty(estim_params_.param_vals)
            indpmodel = estim_params_.param_vals(:,1); %values correspond to parameters declaration order, row number corresponds to order in estimated_params
        end
        indpstderr=[]; %initialize index for stderr parameters
        if ~isempty(estim_params_.var_exo)
            indpstderr = estim_params_.var_exo(:,1); %values correspond to varexo declaration order, row number corresponds to order in estimated_params
        end
        indpcorr=[]; %initialize matrix for corr paramters
        if ~isempty(estim_params_.corrx)
            indpcorr = estim_params_.corrx(:,1:2); %values correspond to varexo declaration order, row number corresponds to order in estimated_params
        end
        if estim_params_.nvn || estim_params_.ncn %nvn is number of stderr parameters and ncn is number of corr parameters of measurement innovations as declared in estimated_params
            error('Analytic computation of standard errrors does not (yet) support measurement errors.\nInstead, define them explicitly as varexo and provide measurement equations in the model definition.\nAlternatively, use numerical standard errors.')
        end
        modparam_nbr = estim_params_.np;        % number of model parameters as declared in estimated_params
        stderrparam_nbr = estim_params_.nvx;    % number of stderr parameters
        corrparam_nbr = estim_params_.ncx;      % number of corr parameters
        totparam_nbr = stderrparam_nbr+corrparam_nbr+modparam_nbr;
        dr.derivs = get_perturbation_params_derivs(M_, options_mom_, estim_params_, oo_, indpmodel, indpstderr, indpcorr, 0); %analytic derivatives of perturbation matrices
        oo_.mom.model_moments_params_derivs = NaN(options_mom_.mom.mom_nbr,totparam_nbr);
        pruned_state_space = pruned_state_space_system(M_, options_mom_, dr, oo_.dr.obs_var, options_mom_.ar, 0, 1);
    else
        pruned_state_space = pruned_state_space_system(M_, options_mom_, dr, oo_.dr.obs_var, options_mom_.ar, 0, 0);
    end
    
    oo_.mom.model_moments = NaN(options_mom_.mom.mom_nbr,1);
    for jm = 1:size(M_.matched_moments,1)
        % First moments
        if ~options_mom_.prefilter && (sum(M_.matched_moments{jm,3}) == 1)
            idx1 = (oo_.dr.obs_var == find(oo_.dr.order_var==M_.matched_moments{jm,1}) );
            oo_.mom.model_moments(jm,1) = pruned_state_space.E_y(idx1);
            if options_mom_.mom.compute_derivs && ( options_mom_.mom.analytic_standard_errors || options_mom_.mom.analytic_jacobian )
                oo_.mom.model_moments_params_derivs(jm,:) = pruned_state_space.dE_y(idx1,:);
            end
        end
        % Second moments
        if (sum(M_.matched_moments{jm,3}) == 2)
            idx1 = (oo_.dr.obs_var == find(oo_.dr.order_var==M_.matched_moments{jm,1}(1)) );
            idx2 = (oo_.dr.obs_var == find(oo_.dr.order_var==M_.matched_moments{jm,1}(2)) );
            if nnz(M_.matched_moments{jm,2}) == 0
                % Covariance
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
                % Autocovariance
                lag = -M_.matched_moments{jm,2}(2); %note that leads/lags in matched_moments are transformed such that first entry is always 0 and the second is a lag
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


elseif strcmp(options_mom_.mom.mom_method,'SMM')
    %------------------------------------------------------------------------------
    % 3. Compute Moments of the model solution for normal innovations
    %------------------------------------------------------------------------------
    
    % create shock series with correct covariance matrix from iid standard normal shocks
    i_exo_var = setdiff(1:M_.exo_nbr, find(diag(M_.Sigma_e) == 0 )); %find singular entries in covariance
    chol_S = chol(M_.Sigma_e(i_exo_var,i_exo_var));
    scaled_shock_series = zeros(size(options_mom_.mom.shock_series)); %initialize
    scaled_shock_series(:,i_exo_var) = options_mom_.mom.shock_series(:,i_exo_var)*chol_S; %set non-zero entries
    
    % simulate series
    y_sim = simult_(M_, options_mom_, dr.ys, dr, scaled_shock_series, options_mom_.order);
    % provide meaningful penalty if data is nan or inf
    if any(any(isnan(y_sim))) || any(any(isinf(y_sim)))
        if options_mom_.vector_output == 1 % lsqnonlin requires vector output
            fval = Inf(size(oo_.mom.Sw,1),1);
        else
            fval = Inf;
        end
        info(1)=180;
        info(4) = 0.1;
        exit_flag = 0;
        if options_mom_.vector_output == 1 % lsqnonlin requires vector output
            fval = ones(size(oo_.mom.data_moments,1),1)*options_mom_.huge_number;
        end
        return
    end
    
    % Remove burn-in and focus on observables (note that y_sim is in declaration order)
    y_sim = y_sim(oo_.dr.order_var(oo_.dr.obs_var) , end-options_mom_.mom.long+1:end)';
    
    if ~all(diag(M_.H)==0)
        i_ME = setdiff([1:size(M_.H,1)],find(diag(M_.H) == 0)); % find ME with 0 variance
        chol_S = chol(M_.H(i_ME,i_ME)); %decompose rest
        shock_mat=zeros(size(options_mom_.mom.ME_shock_series)); %initialize
        shock_mat(:,i_ME)=options_mom_.mom.ME_shock_series(:,i_ME)*chol_S;
        y_sim = y_sim+shock_mat;
    end

    % Remove mean if centered moments
    if options_mom_.prefilter
        y_sim = bsxfun(@minus, y_sim, mean(y_sim,1));
    end
    oo_.mom.model_moments = mom.data_moments(y_sim, oo_, M_.matched_moments, options_mom_);
    
end

%--------------------------------------------------------------------------
% 4. Compute quadratic target function
%--------------------------------------------------------------------------
moments_difference = oo_.mom.data_moments - oo_.mom.model_moments;
residuals = sqrt(options_mom_.mom.weighting_matrix_scaling_factor)*oo_.mom.Sw*moments_difference;
oo_.mom.Q = residuals'*residuals;
if options_mom_.vector_output == 1 % lsqnonlin requires vector output
    fval = residuals;
    if options_mom_.mom.penalized_estimator
        fval=[fval;(xparam1-oo_.prior.mean)./sqrt(diag(oo_.prior.variance))];
    end
else    
    fval = oo_.mom.Q;
    if options_mom_.mom.penalized_estimator
        fval=fval+(xparam1-oo_.prior.mean)'/oo_.prior.variance*(xparam1-oo_.prior.mean);
    end
end

if options_mom_.mom.compute_derivs && options_mom_.mom.analytic_jacobian
    if options_mom_.mom.penalized_estimator
        dxparam1 = eye(length(xparam1));
    end
    
    for jp=1:length(xparam1)
        dmoments_difference = - oo_.mom.model_moments_params_derivs(:,jp);
        dresiduals = sqrt(options_mom_.mom.weighting_matrix_scaling_factor)*oo_.mom.Sw*dmoments_difference;
        
        if options_mom_.vector_output == 1 % lsqnonlin requires vector output            
            if options_mom_.mom.penalized_estimator                
                df(:,jp)=[dresiduals;dxparam1(:,jp)./sqrt(diag(oo_.prior.variance))];
            else
                df(:,jp) = dresiduals;
            end
        else
            df(:,jp) = dresiduals'*residuals + residuals'*dresiduals;
            if options_mom_.mom.penalized_estimator
                df(:,jp)=df(:,jp)+(dxparam1(:,jp))'/oo_.prior.variance*(xparam1-oo_.prior.mean)+(xparam1-oo_.prior.mean)'/oo_.prior.variance*(dxparam1(:,jp));
            end
        end
    end
end


end%main function end

