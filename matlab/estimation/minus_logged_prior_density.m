function [fval, info, exitflag, fake1, fake2] = minus_logged_prior_density(xparams, Prior, options_, M_, estim_params_, oo_)

% Evaluates minus the logged prior density.
 %
% INPUTS
% - xparams            [double]   vector of parameters.
% - Prior              [dprior]   vector specifying prior densities shapes.
% - DynareOptions      [struct]   Options, AKA options_
% - DynareModel        [struct]   Model description, AKA M_
% - EstimatedParams    [struct]   Info about estimated parameters, AKA estimated_params_
% - DynareResults      [struct]   Results, AKA oo_
%
% OUTPUTS
% - fval               [double]  value of minus the logged prior density.
% - info               [double]  4×1 vector, second entry stores penalty, first entry the error code, last entry a penalty (used for optimization).

% Copyright © 2009-2023 Dynare Team
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

exitflag = true;
info = zeros(4,1);
fake1 = [];
fake2 = [];

%------------------------------------------------------------------------------
% 1. Get the structural parameters & define penalties
%------------------------------------------------------------------------------

% Return, with endogenous penalty, if some parameters are smaller than the lower bound of the prior domain.
if ~isequal(options_.mode_compute, 1) && any(xparams<Prior.p3)
    k = find(xparams<Prior.p3);
    fval = Inf;
    exitflag = false;
    info(1) = 41;
    info(4) = sum((Prior.p3(k)-xparams(k)).^2);
    return
end

% Return, with endogenous penalty, if some parameters are greater than the upper bound of the prior domain.
if ~isequal(options_.mode_compute, 1) && any(xparams>Prior.p4)
    k = find(xparams>Prior.p4);
    fval = Inf;
    exitflag = false;
    info(1) = 42;
    info(4) = sum((xparams(k)-Prior.p4(k)).^2);
    return
end

% Get the diagonal elements of the covariance matrices for the structural innovations (Q) and the measurement error (H).
M_ = set_all_parameters(xparams, estim_params_, M_);

Q = M_.Sigma_e;
H = M_.H;

% Test if Q is positive definite.
if ~issquare(Q) || estim_params_.ncx || isfield(estim_params_, 'calibrated_covariances')
    % Try to compute the cholesky decomposition of Q (possible iff Q is positive definite)
    [Q_is_positive_definite, penalty] = ispd(Q);
    if ~Q_is_positive_definite
        % The variance-covariance matrix of the structural innovations is not definite positive. We have to compute the
        % eigenvalues of this matrix in order to build the endogenous penalty.
        fval = Inf;
        exitflag = false;
        info(1) = 43;
        info(4) = penalty;
        return
    end
    if isfield(estim_params_, 'calibrated_covariances')
        correct_flag = check_consistency_covariances(Q);
        if ~correct_flag
            penalty = sum(Q(estim_params_.calibrated_covariances.position).^2);
            fval = Inf;
            exitflag = false;
            info(1) = 71;
            info(4) = penalty;
            return4
        end
    end
end

% Test if H is positive definite.
if ~issquare(H) || estim_params_.ncn || isfield(estim_params_, 'calibrated_covariances_ME')
    [H_is_positive_definite, penalty] = ispd(H);
    if ~H_is_positive_definite
        % The variance-covariance matrix of the measurement errors is not definite positive. We have to compute the eigenvalues
        % of this matrix in order to build the endogenous penalty.
        fval = Inf;
        exitflag = false;
        info(1) = 44;
        info(4) = penalty;
        return
    end
    if isfield(estim_params_, 'calibrated_covariances_ME')
        correct_flag = check_consistency_covariances(H);
        if ~correct_flag
            penalty = sum(H(estim_params_.calibrated_covariances_ME.position).^2);
            fval = Inf;
            exitflag = false;
            info(1) = 72;
            info(4) = penalty;
            return
        end
    end
end


%-----------------------------
% 2. Check BK and steady state
%-----------------------------

[~, info] = resol(0, M_, options_, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);

% Return, with endogenous penalty when possible, if dynare_resolve issues an error code (defined in resol).
if info(1)
    if info(1) == 3 || info(1) == 4 || info(1) == 5 || info(1)==6 ||info(1) == 19 ||...
                info(1) == 20 || info(1) == 21 || info(1) == 23 || info(1) == 26 || ...
                info(1) == 81 || info(1) == 84 ||  info(1) == 85
        %meaningful second entry of output that can be used
        fval = Inf;
        info(4) = info(2);
        exitflag = false;
        return
    else
        fval = Inf;
        info(4) = 0.1;
        exitflag = false;
        return
    end
end

fval = - Prior.density(xparams);
