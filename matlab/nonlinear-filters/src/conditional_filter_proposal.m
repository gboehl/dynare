function [ProposalStateVector, Weights, flag] = conditional_filter_proposal(ReducedForm, y, StateVectors, SampleWeights, Q_lower_triangular_cholesky, H_lower_triangular_cholesky, ...
                                                  H, ParticleOptions, ThreadsOptions, DynareOptions, Model)

% Computes the proposal for each past particle using Gaussian approximations
% for the state errors and the Kalman filter
%
% INPUTS
% - ReducedForm                    [structure]    Matlab's structure describing the reduced form model.
% - y                              [double]       p×1 vector, current observation (p is the number of observed variables).
% - StateVectors
% - SampleWeights
% - Q_lower_triangular_cholesky
% - H_lower_triangular_cholesky
% - H
% - ParticleOptions
% - ThreadsOptions
% - DynareOptions
% - Model
%
% OUTPUTS
% - ProposalStateVector
% - Weights
% - flag

% Copyright © 2012-2020 Dynare Team
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

flag = false;

if ReducedForm.use_k_order_solver
    dr = ReducedForm.dr;
else
    % Set local state space model (first-order approximation).
    ghx  = ReducedForm.ghx;
    ghu  = ReducedForm.ghu;
    % Set local state space model (second-order approximation).
    ghxx = ReducedForm.ghxx;
    ghuu = ReducedForm.ghuu;
    ghxu = ReducedForm.ghxu;
end

constant = ReducedForm.constant;
state_variables_steady_state = ReducedForm.state_variables_steady_state;

mf0 = ReducedForm.mf0;
mf1 = ReducedForm.mf1;
number_of_state_variables = length(mf0);
number_of_observed_variables = length(mf1);
number_of_structural_innovations = length(ReducedForm.Q);

if ParticleOptions.proposal_approximation.montecarlo
    nodes = randn(ParticleOptions.number_of_particles/10, number_of_structural_innovations);
    weights = 1.0/ParticleOptions.number_of_particles;
    weights_c = weights;
elseif ParticleOptions.proposal_approximation.cubature
    [nodes, weights] = spherical_radial_sigma_points(number_of_structural_innovations);
    weights_c = weights;
elseif ParticleOptions.proposal_approximation.unscented
    [nodes, weights, weights_c] = unscented_sigma_points(number_of_structural_innovations, ParticleOptions);
else
    error('Estimation: This approximation for the proposal is not implemented or unknown!')
end

epsilon = Q_lower_triangular_cholesky*nodes';
yhat = repmat(StateVectors-state_variables_steady_state, 1, size(epsilon, 2));

if ReducedForm.use_k_order_solver
    tmp = local_state_space_iteration_k(yhat, epsilon, dr, Model, DynareOptions);
else
    tmp = local_state_space_iteration_2(yhat, epsilon, ghx, ghu, constant, ghxx, ghuu, ghxu, ThreadsOptions.local_state_space_iteration_2);
end

PredictedStateMean = tmp(mf0,:)*weights;
PredictedObservedMean = tmp(mf1,:)*weights;

if ParticleOptions.proposal_approximation.cubature || ParticleOptions.proposal_approximation.montecarlo
    PredictedStateMean = sum(PredictedStateMean, 2);
    PredictedObservedMean = sum(PredictedObservedMean, 2);
    dState = bsxfun(@minus, tmp(mf0,:), PredictedStateMean)'.*sqrt(weights);
    dObserved = bsxfun(@minus, tmp(mf1,:), PredictedObservedMean)'.*sqrt(weights);
    PredictedStateVariance = dState*dState';
    big_mat = [dObserved dState; H_lower_triangular_cholesky zeros(number_of_observed_variables,number_of_state_variables)];
    [~, mat] = qr2(big_mat,0);
    mat = mat';
    PredictedObservedVarianceSquareRoot = mat(1:number_of_observed_variables, 1:number_of_observed_variables);
    CovarianceObservedStateSquareRoot = mat(number_of_observed_variables+(1:number_of_state_variables),1:number_of_observed_variables);
    StateVectorVarianceSquareRoot = mat(number_of_observed_variables+(1:number_of_state_variables),number_of_observed_variables+(1:number_of_state_variables));
    Error = y-PredictedObservedMean;
    StateVectorMean = PredictedStateMean+(CovarianceObservedStateSquareRoot/PredictedObservedVarianceSquareRoot)*Error;
    if ParticleOptions.cpf_weights_method.amisanotristani
        Weights = SampleWeights.*probability2(zeros(number_of_observed_variables,1), PredictedObservedVarianceSquareRoot, Error);
    end
else
    dState = bsxfun(@minus, tmp(mf0,:), PredictedStateMean);
    dObserved = bsxfun(@minus, tmp(mf1,:), PredictedObservedMean);
    PredictedStateVariance = dState*diag(weights_c)*dState';
    PredictedObservedVariance = dObserved*diag(weights_c)*dObserved'+H;
    PredictedStateAndObservedCovariance = dState*diag(weights_c)*dObserved';
    KalmanFilterGain = PredictedStateAndObservedCovariance/PredictedObservedVariance;
    Error = y-PredictedObservedMean;
    StateVectorMean = PredictedStateMean+KalmanFilterGain*Error;
    StateVectorVariance = PredictedStateVariance-KalmanFilterGain*PredictedObservedVariance*KalmanFilterGain';
    StateVectorVariance = 0.5*(StateVectorVariance+StateVectorVariance');
    [StateVectorVarianceSquareRoot, p] = chol(StateVectorVariance, 'lower') ;
    if p
        flag = true;
        ProposalStateVector = zeros(number_of_state_variables, 1);
        Weights = 0.0;
        return
    end
    if ParticleOptions.cpf_weights_method.amisanotristani
        Weights = SampleWeights.*probability2(zeros(number_of_observed_variables, 1), chol(PredictedObservedVariance)', Error);
    end
end

ProposalStateVector = StateVectorVarianceSquareRoot*randn(size(StateVectorVarianceSquareRoot, 2), 1)+StateVectorMean;
if ParticleOptions.cpf_weights_method.murrayjonesparslow
    PredictedStateVariance = 0.5*(PredictedStateVariance+PredictedStateVariance');
    [PredictedStateVarianceSquareRoot, p] = chol(PredictedStateVariance, 'lower');
    if p
        flag = true;
        ProposalStateVector = zeros(number_of_state_variables,1);
        Weights = 0.0;
        return
    end
    Prior = probability2(PredictedStateMean, PredictedStateVarianceSquareRoot, ProposalStateVector);
    Posterior = probability2(StateVectorMean, StateVectorVarianceSquareRoot, ProposalStateVector);
    Likelihood = probability2(y, H_lower_triangular_cholesky, measurement_equations(ProposalStateVector, ReducedForm, ThreadsOptions, DynareOptions, Model));
    Weights = SampleWeights.*Likelihood.*(Prior./Posterior);
end
