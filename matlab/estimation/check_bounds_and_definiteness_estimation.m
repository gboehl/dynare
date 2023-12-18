function [fval,info,exit_flag,Q,H]=check_bounds_and_definiteness_estimation(xparam1, M_, estim_params_, bounds)
% [fval,info,exit_flag,Q,H]=check_bounds_and_definiteness_estimation(xparam1, M_, estim_params_, bounds)
% Checks whether parameter vector satisfies bounds and the positive
% definiteness of covariance matrices
%
% INPUTS
% - xparam1                 [double]              n by 1 vector, estimated parameters.
% - M_                      [struct]              Matlab's structure describing the Model.
% - estim_params_           [struct]              Matlab's structure describing the estimated_parameters.
% - bounds                  [struct]              Matlab's structure specifying the bounds on the paramater values (initialized by dynare_estimation_init).
%
% OUTPUTS
% - fval                    [double]              scalar, value of the likelihood or posterior kernel.
% - info                    [integer]             4 by 1 vector, informations resolution of the model and evaluation of the likelihood.
% - exit_flag               [integer]             scalar, equal to 1 (no issues when evaluating the likelihood) or 0 (not able to evaluate the likelihood).
% - Q                       [matrix]              Covariance matrix of structural shocks
% - H                       [matrix]              Covariance matrix of measurement errors

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

fval        = [];
exit_flag   = 1;
info        = zeros(4,1);
Q=[];
H=[];
% Return, with endogenous penalty, if some parameters are smaller than the lower bound of the prior domain.
if any(xparam1<bounds.lb)
    k = find(xparam1(:) < bounds.lb);
    fval = Inf;
    exit_flag = 0;
    info(1) = 41;
    info(4) = sum((bounds.lb(k)-xparam1(k)).^2);
    return
end

% Return, with endogenous penalty, if some parameters are greater than the upper bound of the prior domain.
if any(xparam1>bounds.ub)
    k = find(xparam1(:)>bounds.ub);
    fval = Inf;
    exit_flag = 0;
    info(1) = 42;
    info(4) = sum((xparam1(k)-bounds.ub(k)).^2);
    return
end

Q = M_.Sigma_e;
H = M_.H;

if ~issquare(Q) || estim_params_.ncx || isfield(estim_params_,'calibrated_covariances')
    [Q_is_positive_definite, penalty] = ispd(Q(estim_params_.Sigma_e_entries_to_check_for_positive_definiteness,estim_params_.Sigma_e_entries_to_check_for_positive_definiteness));
    if ~Q_is_positive_definite
        fval = Inf;
        exit_flag = 0;
        info(1) = 43;
        info(4) = penalty;
        return
    end
    if isfield(estim_params_,'calibrated_covariances')
        correct_flag=check_consistency_covariances(Q);
        if ~correct_flag
            penalty = sum(Q(estim_params_.calibrated_covariances.position).^2);
            fval = Inf;
            exit_flag = 0;
            info(1) = 71;
            info(4) = penalty;
            return
        end
    end

end

if ~issquare(H) || estim_params_.ncn || isfield(estim_params_,'calibrated_covariances_ME')
    [H_is_positive_definite, penalty] = ispd(H(estim_params_.H_entries_to_check_for_positive_definiteness,estim_params_.H_entries_to_check_for_positive_definiteness));
    if ~H_is_positive_definite
        fval = Inf;
        exit_flag = 0;
        info(1) = 44;
        info(4) = penalty;
        return
    end
    if isfield(estim_params_,'calibrated_covariances_ME')
        correct_flag=check_consistency_covariances(H);
        if ~correct_flag
            penalty = sum(H(estim_params_.calibrated_covariances_ME.position).^2);
            fval = Inf;
            exit_flag = 0;
            info(1) = 72;
            info(4) = penalty;
            return
        end
    end
end
