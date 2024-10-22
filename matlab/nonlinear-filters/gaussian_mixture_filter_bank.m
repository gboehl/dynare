function [StateMuPrior,StateSqrtPPrior,StateWeightsPrior,StateMuPost,StateSqrtPPost,StateWeightsPost] =...
    gaussian_mixture_filter_bank(ReducedForm, obs, StateMu, StateSqrtP, StateWeights, ...
                                 StructuralShocksMu, StructuralShocksSqrtP, StructuralShocksWeights, ...
                                 ObservationShocksWeights, H, H_lower_triangular_cholesky, normfactO, ...
                                 ParticleOptions, ThreadsOptions, options_, M_)

% Computes the proposal with a gaussian approximation for importance
% sampling
% This proposal is a gaussian distribution calculated à la Kalman
%
% INPUTS
%    reduced_form_model     [structure] Matlab's structure describing the reduced form model.
%                                       reduced_form_model.measurement.H   [double]   (pp x pp) variance matrix of measurement errors.
%                                       reduced_form_model.state.Q         [double]   (qq x qq) variance matrix of state errors.
%                                       reduced_form_model.state.dr        [structure] output of resol.m.
%    Y                      [double]    pp*smpl matrix of (detrended) data, where pp is the maximum number of observed variables.
%
% OUTPUTS
%    LIK        [double]    scalar, likelihood
%    lik        [double]    vector, density of observations in each period.
%
% REFERENCES
%
% NOTES
%   The vector "lik" is used to evaluate the jacobian of the likelihood.

% Copyright © 2009-2022 Dynare Team
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

order = options_.order;

if ReducedForm.use_k_order_solver
    dr = ReducedForm.dr;
    udr = ReducedForm.udr;
else
    % Set local state space model (first-order approximation).
    ghx  = ReducedForm.ghx;
    ghu  = ReducedForm.ghu;
    % Set local state space model (second-order approximation).
    ghxx = ReducedForm.ghxx;
    ghuu = ReducedForm.ghuu;
    ghxu = ReducedForm.ghxu;
    ghs2 = ReducedForm.ghs2;
    if order == 3
        % Set local state space model (third order approximation).
        ghxxx = ReducedForm.ghxxx;
        ghuuu = ReducedForm.ghuuu;
        ghxxu = ReducedForm.ghxxu;
        ghxuu = ReducedForm.ghxuu;
        ghxss = ReducedForm.ghxss;
        ghuss = ReducedForm.ghuss;
    end
end

constant = ReducedForm.constant;
steadystate = ReducedForm.steadystate;
state_variables_steady_state = ReducedForm.state_variables_steady_state;

mf0 = ReducedForm.mf0;
mf1 = ReducedForm.mf1;
number_of_state_variables = length(mf0);
number_of_observed_variables = length(mf1);
number_of_structural_innovations = length(ReducedForm.Q);

numb = number_of_state_variables+number_of_structural_innovations;

if ParticleOptions.proposal_approximation.cubature
    [nodes3, weights3] = spherical_radial_sigma_points(numb);
    weights_c3 = weights3;
elseif ParticleOptions.proposal_approximation.unscented
    [nodes3, weights3, weights_c3] = unscented_sigma_points(numb, ParticleOptions);
else
    error('This approximation for the proposal is unknown!')
end

epsilon = bsxfun(@plus, StructuralShocksSqrtP*nodes3(:,number_of_state_variables+1:number_of_state_variables+number_of_structural_innovations)', StructuralShocksMu);
StateVectors = bsxfun(@plus, StateSqrtP*nodes3(:,1:number_of_state_variables)', StateMu);
yhat = bsxfun(@minus, StateVectors, state_variables_steady_state);
if ReducedForm.use_k_order_solver
    tmp = local_state_space_iteration_k(yhat, epsilon, dr, M_, options_, udr);
else
    if order == 2
        tmp = local_state_space_iteration_2(yhat, epsilon, ghx, ghu, constant, ghxx, ghuu, ghxu, ThreadsOptions.local_state_space_iteration_2);
    elseif order == 3
        tmp = local_state_space_iteration_3(yhat, epsilon, ghx, ghu, ghxx, ghuu, ghxu, ghs2, ghxxx, ghuuu, ghxxu, ghxuu, ghxss, ghuss, steadystate, ThreadsOptions.local_state_space_iteration_3, false);
    else
        error('Order > 3: use_k_order_solver should be set to true');
    end
end
PredictedStateMean = tmp(mf0,:)*weights3;
PredictedObservedMean = tmp(mf1,:)*weights3;

if ParticleOptions.proposal_approximation.cubature
    PredictedStateMean = sum(PredictedStateMean, 2);
    PredictedObservedMean = sum(PredictedObservedMean, 2);
    dState = (bsxfun(@minus, tmp(mf0,:), PredictedStateMean)').*sqrt(weights3);
    dObserved = (bsxfun(@minus, tmp(mf1,:), PredictedObservedMean)').*sqrt(weights3);
    PredictedStateVariance = dState'*dState;
    big_mat = [dObserved,  dState ; H_lower_triangular_cholesky, zeros(number_of_observed_variables, number_of_state_variables)];
    [~, mat] = qr2(big_mat, 0);
    mat = mat';
    PredictedObservedVarianceSquareRoot = mat(1:number_of_observed_variables, 1:number_of_observed_variables);
    CovarianceObservedStateSquareRoot = mat(number_of_observed_variables+(1:number_of_state_variables), 1:number_of_observed_variables);
    StateVectorVarianceSquareRoot = mat(number_of_observed_variables+(1:number_of_state_variables), number_of_observed_variables+(1:number_of_state_variables));
    iPredictedObservedVarianceSquareRoot = inv(PredictedObservedVarianceSquareRoot);
    iPredictedObservedVariance = iPredictedObservedVarianceSquareRoot'*iPredictedObservedVarianceSquareRoot;
    sqrdet = 1/sqrt(det(iPredictedObservedVariance));
    PredictionError = obs - PredictedObservedMean;
    StateVectorMean = PredictedStateMean + CovarianceObservedStateSquareRoot*iPredictedObservedVarianceSquareRoot*PredictionError;
else
    dState = bsxfun(@minus, tmp(mf0,:), PredictedStateMean);
    dObserved = bsxfun(@minus, tmp(mf1,:), PredictedObservedMean);
    PredictedStateVariance = dState*diag(weights_c3)*dState';
    PredictedObservedVariance = dObserved*diag(weights_c3)*dObserved' + H;
    PredictedStateAndObservedCovariance = dState*diag(weights_c3)*dObserved';
    sqrdet = sqrt(det(PredictedObservedVariance));
    iPredictedObservedVariance = inv(PredictedObservedVariance);
    PredictionError = obs - PredictedObservedMean;
    KalmanFilterGain = PredictedStateAndObservedCovariance*iPredictedObservedVariance;
    StateVectorMean = PredictedStateMean + KalmanFilterGain*PredictionError;
    StateVectorVariance = PredictedStateVariance - KalmanFilterGain*PredictedObservedVariance*KalmanFilterGain';
    StateVectorVariance = .5*(StateVectorVariance+StateVectorVariance');
    StateVectorVarianceSquareRoot = reduced_rank_cholesky(StateVectorVariance)';
end

data_lik_GM_g = exp(-0.5*PredictionError'*iPredictedObservedVariance*PredictionError)/abs(normfactO*sqrdet) + 1e-99;
StateMuPrior = PredictedStateMean;
StateSqrtPPrior = reduced_rank_cholesky(PredictedStateVariance)';
StateWeightsPrior = StateWeights*StructuralShocksWeights;
StateMuPost = StateVectorMean;
StateSqrtPPost = StateVectorVarianceSquareRoot;
StateWeightsPost = StateWeightsPrior*ObservationShocksWeights*data_lik_GM_g;
