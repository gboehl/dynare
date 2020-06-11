function [fval, info, exit_flag, DynareResults, Model, OptionsMoM] = method_of_moments_GMM(xparam1, Bounds, DynareResults, EstimatedParameters, MatchedMoments, Model, OptionsMoM)
% [fval, info, exit_flag, DynareResults, Model, OptionsMoM] = method_of_moments_GMM(xparam1, Bounds, DynareResults, EstimatedParameters, MatchedMoments, Model, OptionsMoM)
% -------------------------------------------------------------------------
% This function evaluates the objective function for GMM estimation
% =========================================================================
% INPUTS
%   o xparam1:                  current value of estimated parameters as returned by set_prior()
%   o Bounds:                   structure containing parameter bounds
%   o DynareResults:            structure for results (oo_)
%   o EstimatedParameters:      structure describing the estimated_parameters (estim_params_)
%   o MatchedMoments:           structure containing information about selected moments to match in estimation (matched_moments_)
%   o Model                     structure describing the Model
%   o OptionsMoM:               structure information about all settings (specified by the user, preprocessor, and taken from global options_)
% -------------------------------------------------------------------------
% OUTPUTS
%   o fval:                     value of the quadratic form of the moment difference (except for lsqnonlin, where this is done implicitly)
%   o info:                     vector storing error code and penalty 
%   o exit_flag:                0 if no error, 1 of error
%   o DynareResults:            structure containing the results with the following updated fields:
%      - mom.modelMoments       [numMom x 1] vector with model moments
%      - mom.Q                  value of the quadratic form of the moment difference
%   o Model:                    Matlab's structure describing the Model
% -------------------------------------------------------------------------
% This function is called by
%  o driver.m
% -------------------------------------------------------------------------
% This function calls
%  o ispd
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
% To Do: check penalized estimation for different optimizers, what is special about mode_compute=1 [@wmutschl]

%--------------------------------------------------------------------------
% 0. Initialization of the returned variables and others...
%--------------------------------------------------------------------------
exit_flag = 1;
info      = zeros(4,1);
xparam1   = xparam1(:); % Ensure that xparam1 is a column vector

%--------------------------------------------------------------------------
% 1. Get the structural parameters & define penalties
%--------------------------------------------------------------------------

% Return, with endogenous penalty, if some parameters are smaller than the lower bound of the parameters.
if any(xparam1<Bounds.lb)    
    if ~isequal(OptionsMoM.mode_compute,1)
        k = find(xparam1<Bounds.lb);
        fval = Inf;
        exit_flag = 0;
        info(1) = 41;
        info(4)= sum((Bounds.lb(k)-xparam1(k)).^2);    
        return
    elseif OptionsMoM.mode_compute == 13
        fval = ones(size(DynareResults.mom.dataMoments,1),1)*OptionsMoM.huge_number;
        return
    end
end

% Return, with endogenous penalty, if some parameters are greater than the upper bound of the parameters.
if any(xparam1>Bounds.ub)
    if ~isequal(OptionsMoM.mode_compute,1)
        k = find(xparam1>Bounds.ub);
        fval = Inf;
        exit_flag = 0;
        info(1) = 42;
        info(4)= sum((xparam1(k)-Bounds.ub(k)).^2);
        return
    elseif OptionsMoM.mode_compute == 13
        fval = ones(size(DynareResults.mom.dataMoments,1),1)*OptionsMoM.huge_number;
        return
    end
end

% Set all parameters
Model = set_all_parameters(xparam1,EstimatedParameters,Model);

% Test if Q is positive definite.
if ~issquare(Model.Sigma_e) || EstimatedParameters.ncx || isfield(EstimatedParameters,'calibrated_covariances')
    [Q_is_positive_definite, penalty] = ispd(Model.Sigma_e(EstimatedParameters.Sigma_e_entries_to_check_for_positive_definiteness,EstimatedParameters.Sigma_e_entries_to_check_for_positive_definiteness));
    if ~Q_is_positive_definite
        if OptionsMoM.mode_compute == 13
            fval = ones(size(DynareResults.mom.dataMoments,1),1)*OptionsMoM.huge_number;
        else
            fval = Inf;
            exit_flag = 0;
            info(1) = 43;
            info(4) = penalty;
        end
        return
    end
    if isfield(EstimatedParameters,'calibrated_covariances')
        correct_flag=check_consistency_covariances(Model.Sigma_e);
        if ~correct_flag
            if OptionsMoM.mode_compute == 13
                fval = ones(size(DynareResults.mom.dataMoments,1),1)*OptionsMoM.huge_number;
            else
                penalty = sum(Model.Sigma_e(EstimatedParameters.calibrated_covariances.position).^2);
                fval = Inf;
                exit_flag = 0;
                info(1) = 71;
                info(4) = penalty;
                if OptionsMoM.mode_compute == 13
                    fval = ones(size(DynareResults.mom.dataMoments,1),1)*OptionsMoM.huge_number;
                end
            end
            return
        end
    end
end

%--------------------------------------------------------------------------
% 2. call resol to compute steady state and model solution
%--------------------------------------------------------------------------

% Compute linear approximation around the deterministic steady state
[dr, info, Model, OptionsMoM, DynareResults] = resol(0, Model, OptionsMoM, DynareResults);

% Return, with endogenous penalty when possible, if resol issues an error code
if info(1)
    if OptionsMoM.mode_compute == 13
        fval = ones(size(DynareResults.mom.dataMoments,1),1)*OptionsMoM.huge_number;
        return
    else
        if info(1) == 3 || info(1) == 4 || info(1) == 5 || info(1)==6 ||info(1) == 19 ||...
                    info(1) == 20 || info(1) == 21 || info(1) == 23 || info(1) == 26 || ...
                    info(1) == 81 || info(1) == 84 ||  info(1) == 85 ||  info(1) == 86
            %meaningful second entry of output that can be used
            fval = Inf;
            info(4) = info(2);
            exit_flag = 0;
            if OptionsMoM.mode_compute == 13
                fval = ones(size(DynareResults.mom.dataMoments,1),1)*OptionsMoM.huge_number;
            end
            return
        else
            fval = Inf;
            info(4) = 0.1;
            exit_flag = 0;
            if OptionsMoM.mode_compute == 13
                fval = ones(size(DynareResults.mom.dataMoments,1),1)*OptionsMoM.huge_number;
            end
            return
        end
    end
end

%--------------------------------------------------------------------------
% 3. Set up pruned state-space system and compute model moments
%--------------------------------------------------------------------------
pruned_state_space = pruned_state_space_system(Model, OptionsMoM, dr, DynareResults.dr.obs_var, OptionsMoM.ar, 0, 0);

DynareResults.mom.modelMoments = nan(OptionsMoM.mom_nbr,1);
offset = 0;
% First moments
if isfield(OptionsMoM.index,'E_y') && nnz(OptionsMoM.index.E_y) > 0 && ~OptionsMoM.prefilter
    E_y = pruned_state_space.E_y;
    E_y_nbr = nnz(OptionsMoM.index.E_y);
    DynareResults.mom.modelMoments(offset+1:E_y_nbr,1) = E_y(OptionsMoM.index.E_y);
    offset = offset + E_y_nbr;
end
% Second moments
if isfield(OptionsMoM.index,'E_yy') && nnz(OptionsMoM.index.E_yy) > 0
    if OptionsMoM.prefilter
        E_yy = pruned_state_space.Var_y;
    else
        E_yy = pruned_state_space.Var_y + pruned_state_space.E_y*pruned_state_space.E_y';
    end    
    E_yy_nbr = nnz(triu(OptionsMoM.index.E_yy));
    DynareResults.mom.modelMoments(offset+(1:E_yy_nbr),1) = E_yy(triu(OptionsMoM.index.E_yy));
    offset = offset + E_yy_nbr;
end

if isfield(OptionsMoM.index,'E_yyt') && nnz(OptionsMoM.index.E_yyt) > 0
    if OptionsMoM.prefilter
        E_yyt = pruned_state_space.Var_yi;
    else
        E_yyt = pruned_state_space.Var_yi + repmat(pruned_state_space.E_y*pruned_state_space.E_y',[1 1 size(pruned_state_space.Var_yi,3)]);
    end
    E_yyt_nbr = nnz(OptionsMoM.index.E_yyt);
    DynareResults.mom.modelMoments(offset+(1:E_yyt_nbr),1) = E_yyt(OptionsMoM.index.E_yyt);
end


%--------------------------------------------------------------------------
% 4. Compute quadratic target function
%--------------------------------------------------------------------------
moments_difference = DynareResults.mom.dataMoments - DynareResults.mom.modelMoments;
residuals = DynareResults.mom.Sw*moments_difference;
DynareResults.mom.Q = residuals'*residuals;
if OptionsMoM.mode_compute == 13 % lsqnonlin
    fval = residuals;
else    
    fval = DynareResults.mom.Q;
    if OptionsMoM.penalized_estimator
        fval=fval+(xparam1-DynareResults.prior.p1)'/diag(DynareResults.prior.p2)*(xparam1-DynareResults.prior.p1);
    end
end


end%main function end

