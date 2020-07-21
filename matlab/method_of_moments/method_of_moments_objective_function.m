function [fval, info, exit_flag, junk1, junk2, oo_, M_, options_mom_] = method_of_moments_objective_function(xparam1, Bounds, oo_, estim_params_, matched_moments_, M_, options_mom_)
% [fval, info, exit_flag, junk1, junk2, oo_, M_, options_mom_] = method_of_moments_objective_function(xparam1, Bounds, oo_, estim_params_, matched_moments_, M_, options_mom_)
% -------------------------------------------------------------------------
% This function evaluates the objective function for GMM/SMM estimation
% =========================================================================
% INPUTS
%   o xparam1:                  current value of estimated parameters as returned by set_prior()
%   o Bounds:                   structure containing parameter bounds
%   o oo_:                      structure for results
%   o estim_params_:            structure describing the estimated_parameters
%   o matched_moments_:         structure containing information about selected moments to match in estimation
%   o M_                        structure describing the model
%   o options_mom_:             structure information about all settings (specified by the user, preprocessor, and taken from global options_)
% -------------------------------------------------------------------------
% OUTPUTS
%   o fval:                     value of the quadratic form of the moment difference (except for lsqnonlin, where this is done implicitly)
%   o info:                     vector storing error code and penalty 
%   o exit_flag:                0 if error, 1 if no error
%   o junk1:                    empty matrix required for optimizer interface
%   o junk2:                    empty matrix required for optimizer interface
%   o oo_:                      structure containing the results with the following updated fields:
%      - mom.model_moments       [numMom x 1] vector with model moments
%      - mom.Q                   value of the quadratic form of the moment difference
%   o M_:                       Matlab's structure describing the model
% -------------------------------------------------------------------------
% This function is called by
%  o method_of_moments.m
% -------------------------------------------------------------------------
% This function calls
%  o check_bounds_and_definiteness_estimation
%  o pruned_state_space_system
%  o resol
%  o set_all_parameters
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

%------------------------------------------------------------------------------
% 0. Initialization of the returned variables and others...
%------------------------------------------------------------------------------

junk1        = [];
junk2        = [];

%--------------------------------------------------------------------------
% 1. Get the structural parameters & define penalties
%--------------------------------------------------------------------------

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
[dr, info, M_, options_mom_, oo_] = resol(0, M_, options_mom_, oo_);

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
    pruned_state_space = pruned_state_space_system(M_, options_mom_, dr, oo_.dr.obs_var, options_mom_.ar, 0, 0);
    
    oo_.mom.model_moments = NaN(options_mom_.mom.mom_nbr,1);
    offset = 0;
    % First moments
    if ~options_mom_.prefilter && isfield(options_mom_.mom.index,'E_y') && nnz(options_mom_.mom.index.E_y) > 0
        E_y = pruned_state_space.E_y;
        E_y_nbr = nnz(options_mom_.mom.index.E_y);
        oo_.mom.model_moments(offset+1:E_y_nbr,1) = E_y(options_mom_.mom.index.E_y);
        offset = offset + E_y_nbr;
    end
    % Second moments
    % Contemporaneous covariance
    if isfield(options_mom_.mom.index,'E_yy') && nnz(options_mom_.mom.index.E_yy) > 0
        if options_mom_.prefilter
            E_yy = pruned_state_space.Var_y;
        else
            E_yy = pruned_state_space.Var_y + pruned_state_space.E_y*pruned_state_space.E_y';
        end
        E_yy_nbr = nnz(tril(options_mom_.mom.index.E_yy));
        oo_.mom.model_moments(offset+(1:E_yy_nbr),1) = E_yy(tril(options_mom_.mom.index.E_yy));
        offset = offset + E_yy_nbr;
    end
    % Lead/lags covariance
    if isfield(options_mom_.mom.index,'E_yyt') && nnz(options_mom_.mom.index.E_yyt) > 0
        if options_mom_.prefilter
            E_yyt = pruned_state_space.Var_yi;
        else
            E_yyt = pruned_state_space.Var_yi + repmat(pruned_state_space.E_y*pruned_state_space.E_y',[1 1 size(pruned_state_space.Var_yi,3)]);
        end
        E_yyt_nbr = nnz(options_mom_.mom.index.E_yyt);
        oo_.mom.model_moments(offset+(1:E_yyt_nbr),1) = E_yyt(options_mom_.mom.index.E_yyt);
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
        if options_mom_.mode_compute==13
            fval = Inf(size(oo_.mom.Sw,1),1);
        else
            fval = Inf;
        end
        info(1)=180;
        info(4) = 0.1;
        exit_flag = 0;
        if options_mom_.mode_compute == 13
            fval = ones(size(oo_.mom.dataMoments,1),1)*options_mom_.huge_number;
        end
        return
    end
    
    % Remove burn-in and focus on observables (note that y_sim is in declaration order)
    y_sim = y_sim(oo_.dr.order_var(oo_.dr.obs_var) , end-options_mom_.mom.long+1:end)';
    
    if ~all(diag(M_.H)==0)
        i_ME = setdiff([1:size(M_.H,1)],find(diag(M_.H) == 0)); % find ME with 0 variance
        chol_S = chol(M_.H(i_ME,i_ME)); %decompose rest
        shock_mat=zeros(size(options_mom_.mom.ME_shock_series)); %initialize
        shock_mat(:,i_ME)=options_mom_.mom.ME_shock_series(:,i_exo_var)*chol_S;
        y_sim = y_sim+shock_mat;
    end

    % Remove mean if centered moments
    if options_mom_.prefilter
        y_sim = bsxfun(@minus, y_sim, mean(y_sim,1));
    end
    oo_.mom.model_moments = method_of_moments_data_moments(y_sim, oo_, matched_moments_, options_mom_);
    
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


end%main function end

