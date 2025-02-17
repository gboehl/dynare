function  IncrementalWeights = gaussian_mixture_densities(obs, StateMuPrior, StateSqrtPPrior, StateWeightsPrior, ...
                                                      StateMuPost, StateSqrtPPost, StateWeightsPost, StateParticles, H, ...
                                                      ReducedForm, ThreadsOptions, options_, M_)

% Elements to calculate the importance sampling ratio
%
% INPUTS
%    reduced_form_model     [structure] Matlab's structure describing the reduced form model.
%                                       reduced_form_model.measurement.H   [double]   (pp x pp) variance matrix of measurement errors.
%                                       reduced_form_model.state.Q         [double]   (qq x qq) variance matrix of state errors.
%                                       reduced_form_model.state.dr        [structure] output of resol.m.
%    Y                      [double]    pp*smpl matrix of (detrended) data, where pp is the maximum number of observed variables.
%    start                  [integer]   scalar, likelihood evaluation starts at 'start'.
%    smolyak_accuracy       [integer]   scalar.
%
% OUTPUTS
%    LIK        [double]    scalar, likelihood
%    lik        [double]    vector, density of observations in each period.
%
% REFERENCES
%
% NOTES
%   The vector "lik" is used to evaluate the jacobian of the likelihood.

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

% Compute the density of particles under the prior distribution
[~, ~, prior] = probability(StateMuPrior, StateSqrtPPrior, StateWeightsPrior, StateParticles);
prior = prior';

% Compute the density of particles under the proposal distribution
[~, ~, proposal] = probability(StateMuPost, StateSqrtPPost, StateWeightsPost, StateParticles);
proposal = proposal';

% Compute the density of the current observation conditionally to each particle
yt_t_1_i = measurement_equations(StateParticles, ReducedForm, ThreadsOptions, options_, M_);

% likelihood
likelihood = probability2(obs, sqrt(H), yt_t_1_i);
IncrementalWeights = likelihood.*prior./proposal;
