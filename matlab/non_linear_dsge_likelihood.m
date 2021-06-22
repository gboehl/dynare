function [fval,info,exit_flag,DLIK,Hess,ys,trend_coeff,Model,DynareOptions,BayesInfo,DynareResults] = non_linear_dsge_likelihood(xparam1,DynareDataset,DatasetInfo,DynareOptions,Model,EstimatedParameters,BayesInfo,BoundsInfo,DynareResults)

% Evaluates the posterior kernel of a dsge model using a non linear filter.
%
% INPUTS
% - xparam1                 [double]              n×1 vector, estimated parameters.
% - DynareDataset           [struct]              Matlab's structure containing the dataset (initialized by dynare, aka dataset_).
% - DatasetInfo             [struct]              Matlab's structure describing the dataset (initialized by dynare, aka dataset_info).
% - DynareOptions           [struct]              Matlab's structure describing the options (initialized by dynare, aka options_).
% - Model                   [struct]              Matlab's structure describing the Model (initialized by dynare, aka M_).
% - EstimatedParameters     [struct]              Matlab's structure describing the estimated_parameters (initialized by dynare, aka estim_params_).
% - BayesInfo               [struct]              Matlab's structure describing the priors (initialized by dynare,aka bayesopt_).
% - BoundsInfo              [struct]              Matlab's structure specifying the bounds on the paramater values (initialized by dynare,aka bayesopt_).
% - DynareResults           [struct]              Matlab's structure gathering the results (initialized by dynare,aka oo_).
%
% OUTPUTS
% - fval                    [double]              scalar, value of the likelihood or posterior kernel.
% - info                    [integer]             4×1 vector, informations resolution of the model and evaluation of the likelihood.
% - exit_flag               [integer]             scalar, equal to 1 (no issues when evaluating the likelihood) or 0 (not able to evaluate the likelihood).
% - DLIK                    [double]              Empty array.
% - Hess                    [double]              Empty array.
% - ys                      [double]              Empty array.
% - trend_coeff             [double]              Empty array.
% - Model                   [struct]              Updated Model structure described in INPUTS section.
% - DynareOptions           [struct]              Updated DynareOptions structure described in INPUTS section.
% - BayesInfo               [struct]              See INPUTS section.
% - DynareResults           [struct]              Updated DynareResults structure described in INPUTS section.

% Copyright (C) 2010-2021 Dynare Team
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

% Initialization of the returned arguments.
fval            = [];
ys              = [];
trend_coeff     = [];
exit_flag       = 1;
DLIK            = [];
Hess            = [];

% Ensure that xparam1 is a column vector.
% (Don't do the transformation if xparam1 is empty, otherwise it would become a
%  0×1 matrix, which create issues with older MATLABs when comparing with [] in
%  check_bounds_and_definiteness_estimation)
if ~isempty(xparam1)
    xparam1 = xparam1(:);
end

% Issue an error if loglinear option is used.
if DynareOptions.loglinear
    error('non_linear_dsge_likelihood: It is not possible to use a non linear filter with the option loglinear!')
end

%------------------------------------------------------------------------------
% 1. Get the structural parameters & define penalties
%------------------------------------------------------------------------------

Model = set_all_parameters(xparam1,EstimatedParameters,Model);

[fval,info,exit_flag,Q,H]=check_bounds_and_definiteness_estimation(xparam1, Model, EstimatedParameters, BoundsInfo);
if info(1)
    return
end

%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------

% Linearize the model around the deterministic sdteadystate and extract the matrices of the state equation (T and R).
[dr, info, Model, DynareResults] = resol(0, Model, DynareOptions, DynareResults);

if info(1)
    if info(1) == 3 || info(1) == 4 || info(1) == 5 || info(1)==6 ||info(1) == 19 || ...
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

% Define a vector of indices for the observed variables. Is this really usefull?...
BayesInfo.mf = BayesInfo.mf1;

% Get needed informations for kalman filter routines.
start = DynareOptions.presample+1;
Y = transpose(DynareDataset.data);

%------------------------------------------------------------------------------
% 3. Initial condition of the Kalman filter
%------------------------------------------------------------------------------

mf0 = BayesInfo.mf0;
mf1 = BayesInfo.mf1;
restrict_variables_idx = dr.restrict_var_list;
state_variables_idx = restrict_variables_idx(mf0);
number_of_state_variables = length(mf0);

ReducedForm.steadystate = dr.ys(dr.order_var(restrict_variables_idx));
ReducedForm.constant = ReducedForm.steadystate + .5*dr.ghs2(restrict_variables_idx);
ReducedForm.state_variables_steady_state = dr.ys(dr.order_var(state_variables_idx));
ReducedForm.Q = Q;
ReducedForm.H = H;
ReducedForm.mf0 = mf0;
ReducedForm.mf1 = mf1;

if DynareOptions.k_order_solver && ~(DynareOptions.particle.pruning && DynareOptions.order==2)
    ReducedForm.use_k_order_solver = true;
    ReducedForm.dr = dr;
else
    ReducedForm.use_k_order_solver = false;
    ReducedForm.ghx  = dr.ghx(restrict_variables_idx,:);
    ReducedForm.ghu  = dr.ghu(restrict_variables_idx,:);
    ReducedForm.ghxx = dr.ghxx(restrict_variables_idx,:);
    ReducedForm.ghuu = dr.ghuu(restrict_variables_idx,:);
    ReducedForm.ghxu = dr.ghxu(restrict_variables_idx,:);
end

% Set initial condition.
switch DynareOptions.particle.initialization
  case 1% Initial state vector covariance is the ergodic variance associated to the first order Taylor-approximation of the model.
    StateVectorMean = ReducedForm.constant(mf0);
    [A,B] = kalman_transition_matrix(dr,dr.restrict_var_list,dr.restrict_columns,Model.exo_nbr);
    StateVectorVariance = lyapunov_symm(A, B*Q*B', DynareOptions.lyapunov_fixed_point_tol, ...
                                        DynareOptions.qz_criterium, DynareOptions.lyapunov_complex_threshold, [], DynareOptions.debug);
    StateVectorVariance = StateVectorVariance(mf0,mf0);
  case 2% Initial state vector covariance is a monte-carlo based estimate of the ergodic variance (consistent with a k-order Taylor-approximation of the model).
    StateVectorMean = ReducedForm.constant(mf0);
    old_DynareOptionsperiods = DynareOptions.periods;
    DynareOptions.periods = 5000;
    old_DynareOptionspruning =  DynareOptions.pruning;
    DynareOptions.pruning = DynareOptions.particle.pruning;
    y_ = simult(DynareResults.steady_state, dr,Model,DynareOptions,DynareResults);
    y_ = y_(dr.order_var(state_variables_idx),2001:5000); %state_variables_idx is in dr-order while simult_ is in declaration order
    if any(any(isnan(y_))) ||  any(any(isinf(y_))) && ~ DynareOptions.pruning
        fval = Inf;
        info(1) = 202;
        info(4) = 0.1;
        exit_flag = 0;
        return;        
    end
    StateVectorVariance = cov(y_');       
    DynareOptions.periods = old_DynareOptionsperiods;
    DynareOptions.pruning = old_DynareOptionspruning;
    clear('old_DynareOptionsperiods','y_');
  case 3% Initial state vector covariance is a diagonal matrix (to be used
        % if model has stochastic trends).
    StateVectorMean = ReducedForm.constant(mf0);
    StateVectorVariance = DynareOptions.particle.initial_state_prior_std*eye(number_of_state_variables);
  otherwise
    error('Unknown initialization option!')
end
ReducedForm.StateVectorMean = StateVectorMean;
ReducedForm.StateVectorVariance = StateVectorVariance;

[~, flag] = chol(ReducedForm.StateVectorVariance);%reduced_rank_cholesky(ReducedForm.StateVectorVariance)';
if flag
    fval = Inf;
    info(1) = 201;
    info(4) = 0.1;
    exit_flag = 0;    
    return;
end
%------------------------------------------------------------------------------
% 4. Likelihood evaluation
%------------------------------------------------------------------------------
DynareOptions.warning_for_steadystate = 0;
[s1,s2] = get_dynare_random_generator_state();
LIK = feval(DynareOptions.particle.algorithm, ReducedForm, Y, start, DynareOptions.particle, DynareOptions.threads, DynareOptions, Model);
set_dynare_random_generator_state(s1,s2);
if imag(LIK)
    fval = Inf;
    info(1) = 46;
    info(4) = 0.1;
    exit_flag = 0;
    return
elseif isnan(LIK)
    fval = Inf;
    info(1) = 45;
    info(4) = 0.1;
    exit_flag = 0;
    return
else
    likelihood = LIK;
end
DynareOptions.warning_for_steadystate = 1;
% ------------------------------------------------------------------------------
% Adds prior if necessary
% ------------------------------------------------------------------------------
lnprior = priordens(xparam1(:),BayesInfo.pshape,BayesInfo.p6,BayesInfo.p7,BayesInfo.p3,BayesInfo.p4);
fval = (likelihood-lnprior);

if isnan(fval)
    fval = Inf;
    info(1) = 47;
    info(4) = 0.1;
    exit_flag = 0;
    return
end

if ~isreal(fval)
    fval = Inf;
    info(1) = 48;
    info(4) = 0.1;
    exit_flag = 0;
    return
end

if isinf(LIK)
    fval = Inf;
    info(1) = 50;
    info(4) = 0.1;
    exit_flag = 0;
    return
end
