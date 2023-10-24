function [fval,info,exit_flag,fake_1,fake_2] = minus_logged_prior_density(xparams,pshape,p6,p7,p3,p4,options_,M_,estim_params_,oo_)
% [fval,info,exit_flag,fake_1,fake_2] = minus_logged_prior_density(xparams,pshape,p6,p7,p3,p4,options_,M_,estim_params_,oo_)
% Evaluates minus the logged prior density.
%
% INPUTS
%   xparams             [double]   vector of parameters.
%   pshape              [integer]  vector specifying prior densities shapes.
%   p6                  [double]   vector, first hyperparameter.
%   p7                  [double]   vector, second hyperparameter.
%   p3                  [double]   vector, prior's lower bound.
%   p4                  [double]   vector, prior's upper bound.
%   prior_sup_bound     [double]   vector, prior's upper bound.
%   options_            [structure] describing the options
%   M_                  [structure] describing the model
%   estim_params_       [structure] characterizing parameters to be estimated
%   oo_                 [structure] storing the results
%
% OUTPUTS
%   f          [double]  value of minus the logged prior density.
%   info       [double]  vector: second entry stores penalty, first entry the error code.
%
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

fake_1 = [];
fake_2 = [];

exit_flag = 1;
info = zeros(4,1);

%------------------------------------------------------------------------------
% 1. Get the structural parameters & define penalties
%------------------------------------------------------------------------------

% Return, with endogenous penalty, if some parameters are smaller than the lower bound of the prior domain.
if ~isequal(options_.mode_compute,1) && any(xparams<p3)
    k = find(xparams<p3);
    fval = Inf;
    exit_flag = 0;
    info(1) = 41;
    info(4) = sum((p3(k)-xparams(k)).^2);
    return
end

% Return, with endogenous penalty, if some parameters are greater than the upper bound of the prior domain.
if ~isequal(options_.mode_compute,1) && any(xparams>p4)
    k = find(xparams>p4);
    fval = Inf;
    exit_flag = 0;
    info(1) = 42;
    info(4) = sum((xparams(k)-p4(k)).^2);
    return
end

% Get the diagonal elements of the covariance matrices for the structural innovations (Q) and the measurement error (H).
M_ = set_all_parameters(xparams,estim_params_,M_);

Q = M_.Sigma_e;
H = M_.H;

% Test if Q is positive definite.
if ~issquare(Q) || estim_params_.ncx || isfield(estim_params_,'calibrated_covariances')
    % Try to compute the cholesky decomposition of Q (possible iff Q is positive definite)
    [Q_is_positive_definite, penalty] = ispd(Q);
    if ~Q_is_positive_definite
        % The variance-covariance matrix of the structural innovations is not definite positive. We have to compute the eigenvalues of this matrix in order to build the endogenous penalty.
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
            return4
        end
    end

end

% Test if H is positive definite.
if ~issquare(H) || estim_params_.ncn || isfield(estim_params_,'calibrated_covariances_ME')
    [H_is_positive_definite, penalty] = ispd(H);
    if ~H_is_positive_definite
        % The variance-covariance matrix of the measurement errors is not definite positive. We have to compute the eigenvalues of this matrix in order to build the endogenous penalty.
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


%-----------------------------
% 2. Check BK and steady state
%-----------------------------

[~,info] = resol(0,M_,options_,oo_.dr,oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);

% Return, with endogenous penalty when possible, if dynare_resolve issues an error code (defined in resol).
if info(1)
    if info(1) == 3 || info(1) == 4 || info(1) == 5 || info(1)==6 ||info(1) == 19 ||...
                info(1) == 20 || info(1) == 21 || info(1) == 23 || info(1) == 26 || ...
                info(1) == 81 || info(1) == 84 ||  info(1) == 85
        %meaningful second entry of output that can be used
        fval = Inf;
        info(4) = info(2);
        exit_flag = 0;
        return
    else
        fval = Inf;
        info(4) = 0.1;
        exit_flag = 0;
        return
    end
end

fval = - priordens(xparams,pshape,p6,p7,p3,p4);