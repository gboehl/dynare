function [LIK,lik] = conditional_particle_filter(ReducedForm, Y, s, ParticleOptions, ThreadsOptions, DynareOptions, Model)

% Evaluates the likelihood of a non-linear model with a particle filter
%
% INPUTS
% - ReducedForm        [structure]    Matlab's structure describing the reduced form model.
% - Y                  [double]       p×T matrix of (detrended) data, where p is the number of observed variables.
% - s                  [integer]      scalar, likelihood evaluation starts at s (has to be smaller than T, the sample length provided in Y).
% - ParticlesOptions   [struct]
% - ThreadsOptions     [struct]
% - DynareOptions      [struct]
% - Model              [struct]
%
% OUTPUTS
% - LIK                [double]    scalar, likelihood
% - lik                [double]    (T-s+1)×1 vector, density of observations in each period.
%
% REMARKS
% - The proposal is built using the Kalman updating step for each particle.
% - we need draws in the errors distributions
% Whether we use Monte-Carlo draws from a multivariate gaussian distribution
% as in Amisano & Tristani (JEDC 2010).
% Whether we use multidimensional Gaussian sparse grids approximations:
% - a univariate Kronrod-Paterson Gaussian quadrature combined by the Smolyak
% operator (ref: Winschel & Kratzig, 2010).
% - a spherical-radial cubature (ref: Arasaratnam & Haykin, 2009a,2009b).
% - a scaled unscented transform cubature (ref: Julier & Uhlmann 1997, van der
% Merwe & Wan 2003).
%
% Pros:
% - Allows using current observable information in the proposal
% - The use of sparse grids Gaussian approximation is much faster than the Monte-Carlo approach
% Cons:
% - The use of the Kalman updating step may biais the proposal distribution since
% it has been derived in a linear context and is implemented in a nonlinear
% context. That is why particle resampling is performed.

% Copyright (C) 2009-2020 Dynare Team
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

% Set default for third input argument.
if isempty(s)
    s = 1;
end

T = size(Y,2);
p = length(ReducedForm.mf1);
n = ParticleOptions.number_of_particles;

% Get covariance matrices
Q = ReducedForm.Q;
H = ReducedForm.H;
if isempty(H)
    H = 0;
    H_lower_triangular_cholesky = 0;
else
    H_lower_triangular_cholesky = chol(H)';
end

% Get initial condition for the state vector.
StateVectorMean = ReducedForm.StateVectorMean;
StateVectorVarianceSquareRoot = chol(ReducedForm.StateVectorVariance)';
state_variance_rank = size(StateVectorVarianceSquareRoot, 2);
Q_lower_triangular_cholesky = chol(Q)';

% Set seed for randn().
set_dynare_seed('default');

% Initialization of the likelihood.
lik = NaN(T, 1);
ks = 0;
StateParticles = bsxfun(@plus, StateVectorVarianceSquareRoot*randn(state_variance_rank, n), StateVectorMean);
SampleWeights = ones(1, n)/n;

for t=1:T
    flags = false(n, 1);
    for i=1:n
        [StateParticles(:,i), SampleWeights(i), flags(i)] = ...
            conditional_filter_proposal(ReducedForm, Y(:,t), StateParticles(:,i), SampleWeights(i), Q_lower_triangular_cholesky, H_lower_triangular_cholesky, H, ParticleOptions, ThreadsOptions, DynareOptions, Model);
    end
    if any(flags)
        LIK = -Inf;
        lik(t) = -Inf;
        return 
    end
    SumSampleWeights = sum(SampleWeights);
    lik(t) = log(SumSampleWeights);
    SampleWeights = SampleWeights./SumSampleWeights;
    if (ParticleOptions.resampling.status.generic && neff(SampleWeights)<ParticleOptions.resampling.threshold*T) || ParticleOptions.resampling.status.systematic
        ks = ks + 1;
        StateParticles = resample(StateParticles', SampleWeights', ParticleOptions)';
        SampleWeights = ones(1, n)/n;
    end
end

LIK = -sum(lik(s:end));