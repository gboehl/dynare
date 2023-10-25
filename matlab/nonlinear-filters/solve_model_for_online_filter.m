function [info, M_, options_, oo_, ReducedForm] = ...
    solve_model_for_online_filter(setinitialcondition, xparam1, dataset_, options_, M_, estim_params_, bayestopt_, bounds, oo_)

% Solves the dsge model for an particular parameters set.
%
% INPUTS
% - setinitialcondition      [logical]    return initial condition if true.
% - xparam1                  [double]     n×1 vector, parameter values.
% - dataset_                 [struct]     Dataset for estimation.
% - options_                 [struct]     Dynare options.
% - M_                       [struct]     Model description.
% - estim_params_            [struct]     Estimated parameters.
% - bayestopt_               [struct]     Prior definition.
% - oo_                      [struct]     Dynare results.
%
% OUTPUTS
% - info                     [integer]    scalar, nonzero if any problem occur when computing the reduced form.
% - M_                       [struct]     M_ description.
% - options_                 [struct]     Dynare options.
% - oo_                      [struct]     Dynare results.
% - ReducedForm              [struct]     Reduced form model.

% Copyright © 2013-2023 Dynare Team
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

persistent init_flag restrict_variables_idx state_variables_idx mf0 mf1 number_of_state_variables

info = 0;

%----------------------------------------------------
% 1. Get the structural parameters & define penalties
%----------------------------------------------------

% Test if some parameters are smaller than the lower bound of the prior domain.
if any(xparam1<bounds.lb)
    info = 41;
    return
end

% Test if some parameters are greater than the upper bound of the prior domain.
if any(xparam1>bounds.ub)
    info = 42;
    return
end

% Get the diagonal elements of the covariance matrices for the structural innovations (Q) and the measurement error (H).
Q = M_.Sigma_e;
H = M_.H;
for i=1:estim_params_.nvx
    k =estim_params_.var_exo(i,1);
    Q(k,k) = xparam1(i)*xparam1(i);
end
offset = estim_params_.nvx;
if estim_params_.nvn
    for i=1:estim_params_.nvn
        H(i,i) = xparam1(i+offset)*xparam1(i+offset);
    end
    offset = offset+estim_params_.nvn;
else
    H = zeros(size(dataset_.data, 2));
end

% Get the off-diagonal elements of the covariance matrix for the structural innovations. Test if Q is positive definite.
if estim_params_.ncx
    for i=1:estim_params_.ncx
        k1 =estim_params_.corrx(i,1);
        k2 =estim_params_.corrx(i,2);
        Q(k1,k2) = xparam1(i+offset)*sqrt(Q(k1,k1)*Q(k2,k2));
        Q(k2,k1) = Q(k1,k2);
    end
    % Try to compute the cholesky decomposition of Q (possible iff Q is positive definite)
    [~, testQ] = chol(Q);
    if testQ
        % The variance-covariance matrix of the structural innovations is not definite positive.
        info = 43;
        return
    end
    offset = offset+estim_params_.ncx;
end

% Get the off-diagonal elements of the covariance matrix for the measurement errors. Test if H is positive definite.
if estim_params_.ncn
    corrn_observable_correspondence = estim_params_.corrn_observable_correspondence;
    for i=1:estim_params_.ncn
        k1 = corrn_observable_correspondence(i,1);
        k2 = corrn_observable_correspondence(i,2);
        H(k1,k2) = xparam1(i+offset)*sqrt(H(k1,k1)*H(k2,k2));
        H(k2,k1) = H(k1,k2);
    end
    % Try to compute the cholesky decomposition of H (possible iff H is positive definite)
    [~, testH] = chol(H);
    if testH
        % The variance-covariance matrix of the measurement errors is not definite positive.
        info = 44;
        return
    end
    offset = offset+estim_params_.ncn;
end

% Update estimated structural parameters in Mode.params.
if estim_params_.np > 0
    M_.params(estim_params_.param_vals(:,1)) = xparam1(offset+1:end);
end

% Update M_.Sigma_e and M_.H.
M_.Sigma_e = Q;
M_.H = H;

%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------

warning('off', 'MATLAB:nearlySingularMatrix')
[~, ~, ~, info, oo_.dr, M_.params] = ...
    dynare_resolve(M_, options_, oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state, 'restrict');
warning('on', 'MATLAB:nearlySingularMatrix')

if info(1)~=0
    if nargout==5
        ReducedForm = 0;
    end
    return
end

% Get decision rules and transition equations.
dr = oo_.dr;

% Set persistent variables (first call).
if isempty(init_flag)
    mf0 = bayestopt_.mf0;
    mf1 = bayestopt_.mf1;
    restrict_variables_idx  = dr.restrict_var_list;
    state_variables_idx = restrict_variables_idx(mf0);
    number_of_state_variables = length(mf0);
    init_flag = true;
end


% Return reduced form model.
if nargout>4
    ReducedForm.ghx = dr.ghx(restrict_variables_idx,:);
    ReducedForm.ghu = dr.ghu(restrict_variables_idx,:);
    ReducedForm.steadystate = dr.ys(dr.order_var(restrict_variables_idx));
    if options_.order==2
        ReducedForm.use_k_order_solver = false;
        ReducedForm.ghxx = dr.ghxx(restrict_variables_idx,:);
        ReducedForm.ghuu = dr.ghuu(restrict_variables_idx,:);
        ReducedForm.ghxu = dr.ghxu(restrict_variables_idx,:);
        ReducedForm.constant = ReducedForm.steadystate + .5*dr.ghs2(restrict_variables_idx);
        ReducedForm.ghs2 = dr.ghs2(restrict_variables_idx,:);
    elseif options_.order>=3
        ReducedForm.use_k_order_solver = true;
        ReducedForm.dr = dr;
        ReducedForm.udr = folded_to_unfolded_dr(dr, M_, options_);
    else
        n_states=size(dr.ghx,2);
        n_shocks=size(dr.ghu,2);
        ReducedForm.use_k_order_solver = false;
        ReducedForm.ghxx = zeros(size(restrict_variables_idx,1),n_states^2);
        ReducedForm.ghuu = zeros(size(restrict_variables_idx,1),n_shocks^2);
        ReducedForm.ghxu = zeros(size(restrict_variables_idx,1),n_states*n_shocks);
        ReducedForm.constant = ReducedForm.steadystate;
    end
    ReducedForm.state_variables_steady_state = dr.ys(dr.order_var(state_variables_idx));
    ReducedForm.Q = Q;
    ReducedForm.H = H;
    ReducedForm.mf0 = mf0;
    ReducedForm.mf1 = mf1;
end

% Set initial condition
if setinitialcondition
    switch options_.particle.initialization
      case 1% Initial state vector covariance is the ergodic variance associated to the first order Taylor-approximation of the model.
        StateVectorMean = ReducedForm.state_variables_steady_state;%.constant(mf0);
        [A,B] = kalman_transition_matrix(dr,dr.restrict_var_list,dr.restrict_columns);
        StateVectorVariance = lyapunov_symm(A, B*ReducedForm.Q*B', options_.lyapunov_fixed_point_tol, ...
                                        options_.qz_criterium, options_.lyapunov_complex_threshold, [], options_.debug);
        StateVectorVariance = StateVectorVariance(mf0,mf0);
      case 2% Initial state vector covariance is a monte-carlo based estimate of the ergodic variance (consistent with a k-order Taylor-approximation of the model).
        StateVectorMean = ReducedForm.state_variables_steady_state;%.constant(mf0);
        old_DynareOptionsperiods = options_.periods;
        options_.periods = 5000;
        old_DynareOptionspruning =  options_.pruning;
        options_.pruning = options_.particle.pruning;
        y_ = simult(dr.ys, dr, M_, options_);
        y_ = y_(dr.order_var(state_variables_idx),2001:options_.periods);
        StateVectorVariance = cov(y_');
        options_.periods = old_DynareOptionsperiods;
        options_.pruning = old_DynareOptionspruning;
        clear('old_DynareOptionsperiods','y_');
      case 3% Initial state vector covariance is a diagonal matrix.
        StateVectorMean = ReducedForm.state_variables_steady_state;%.constant(mf0);
        StateVectorVariance = options_.particle.initial_state_prior_std*eye(number_of_state_variables);
      otherwise
        error('Unknown initialization option!')
    end
    ReducedForm.StateVectorMean = StateVectorMean;
    ReducedForm.StateVectorVariance = StateVectorVariance;
end
