function [info, Model, DynareOptions, DynareResults, ReducedForm] = ...
    solve_model_for_online_filter(setinitialcondition, xparam1, DynareDataset, DynareOptions, Model, EstimatedParameters, BayesInfo, bounds, DynareResults)

% Solves the dsge model for an particular parameters set.
%
% INPUTS
% - setinitialcondition      [logical]    return initial condition if true.
% - xparam1                  [double]     n×1 vector, parameter values.
% - DynareDataset            [struct]     Dataset for estimation (dataset_).
% - DynareOptions            [struct]     Dynare options (options_).
% - Model                    [struct]     Model description (M_).
% - EstimatedParameters      [struct]     Estimated parameters (estim_params_).
% - BayesInfo                [struct]     Prior definition (bayestopt_).
% - DynareResults            [struct]     Dynare results (oo_).
%
% OUTPUTS
% - info                     [integer]    scalar, nonzero if any problem occur when computing the reduced form.
% - Model                    [struct]     Model description (M_).
% - DynareOptions            [struct]     Dynare options (options_).
% - DynareResults            [struct]     Dynare results (oo_).
% - ReducedForm              [struct]     Reduced form model.

% Copyright (C) 2013-2019 Dynare Team
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
Q = Model.Sigma_e;
H = Model.H;
for i=1:EstimatedParameters.nvx
    k =EstimatedParameters.var_exo(i,1);
    Q(k,k) = xparam1(i)*xparam1(i);
end
offset = EstimatedParameters.nvx;
if EstimatedParameters.nvn
    for i=1:EstimatedParameters.nvn
        H(i,i) = xparam1(i+offset)*xparam1(i+offset);
    end
    offset = offset+EstimatedParameters.nvn;
else
    H = zeros(size(DynareDataset.data, 2));
end

% Get the off-diagonal elements of the covariance matrix for the structural innovations. Test if Q is positive definite.
if EstimatedParameters.ncx
    for i=1:EstimatedParameters.ncx
        k1 =EstimatedParameters.corrx(i,1);
        k2 =EstimatedParameters.corrx(i,2);
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
    offset = offset+EstimatedParameters.ncx;
end

% Get the off-diagonal elements of the covariance matrix for the measurement errors. Test if H is positive definite.
if EstimatedParameters.ncn
    corrn_observable_correspondence = EstimatedParameters.corrn_observable_correspondence;
    for i=1:EstimatedParameters.ncn
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
    offset = offset+EstimatedParameters.ncn;
end

% Update estimated structural parameters in Mode.params.
if EstimatedParameters.np > 0
    Model.params(EstimatedParameters.param_vals(:,1)) = xparam1(offset+1:end);
end

% Update Model.Sigma_e and Model.H.
Model.Sigma_e = Q;
Model.H = H;

%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------

warning('off', 'MATLAB:nearlySingularMatrix')
[~, ~, ~, info, Model, DynareResults] = ...
    dynare_resolve(Model, DynareOptions, DynareResults, 'restrict');
warning('on', 'MATLAB:nearlySingularMatrix')

if info(1)~=0
    if nargout==5
        ReducedForm = 0;
    end
    return
end

% Get decision rules and transition equations.
dr = DynareResults.dr;

% Set persistent variables (first call).
if isempty(init_flag)
    mf0 = BayesInfo.mf0;
    mf1 = BayesInfo.mf1;
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
    if DynareOptions.order==2
        ReducedForm.use_k_order_solver = false;
        ReducedForm.ghxx = dr.ghxx(restrict_variables_idx,:);
        ReducedForm.ghuu = dr.ghuu(restrict_variables_idx,:);
        ReducedForm.ghxu = dr.ghxu(restrict_variables_idx,:);
        ReducedForm.constant = ReducedForm.steadystate + .5*dr.ghs2(restrict_variables_idx);
    elseif DynareOptions.order>=3
        ReducedForm.use_k_order_solver = true;
        ReducedForm.dr = dr;
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
    switch DynareOptions.particle.initialization
      case 1% Initial state vector covariance is the ergodic variance associated to the first order Taylor-approximation of the model.
        StateVectorMean = ReducedForm.state_variables_steady_state;%.constant(mf0);
        [A,B] = kalman_transition_matrix(dr,dr.restrict_var_list,dr.restrict_columns,Model.exo_nbr);
        StateVectorVariance2 = lyapunov_symm(ReducedForm.ghx(mf0,:),ReducedForm.ghu(mf0,:)*ReducedForm.Q*ReducedForm.ghu(mf0,:)',DynareOptions.lyapunov_fixed_point_tol,DynareOptions.qz_criterium,DynareOptions.lyapunov_complex_threshold);
        StateVectorVariance = lyapunov_symm(A, B*ReducedForm.Q*B', DynareOptions.lyapunov_fixed_point_tol, ...
                                        DynareOptions.qz_criterium, DynareOptions.lyapunov_complex_threshold, [], DynareOptions.debug);
        StateVectorVariance = StateVectorVariance(mf0,mf0);
      case 2% Initial state vector covariance is a monte-carlo based estimate of the ergodic variance (consistent with a k-order Taylor-approximation of the model).
        StateVectorMean = ReducedForm.state_variables_steady_state;%.constant(mf0);
        old_DynareOptionsperiods = DynareOptions.periods;
        DynareOptions.periods = 5000;
        old_DynareOptionspruning =  DynareOptions.pruning;
        DynareOptions.pruning = DynareOptions.particle.pruning;
        y_ = simult(dr.ys, dr, Model, DynareOptions, DynareResults);
        y_ = y_(dr.order_var(state_variables_idx),2001:DynareOptions.periods);
        StateVectorVariance = cov(y_');
        DynareOptions.periods = old_DynareOptionsperiods;
        DynareOptions.pruning = old_DynareOptionspruning;
        clear('old_DynareOptionsperiods','y_');
      case 3% Initial state vector covariance is a diagonal matrix.
        StateVectorMean = ReducedForm.state_variables_steady_state;%.constant(mf0);
        StateVectorVariance = DynareOptions.particle.initial_state_prior_std*eye(number_of_state_variables);
      otherwise
        error('Unknown initialization option!')
    end
    ReducedForm.StateVectorMean = StateVectorMean;
    ReducedForm.StateVectorVariance = StateVectorVariance;
end