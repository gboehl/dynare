function [forcs, e] = get_shock_paths(cL, H, mcValue, shocks, forcs, T, R, mv, mu)
% [forcs, e] = get_shock_paths(cL, H, mcValue, shocks, forcs, T, R, mv, mu)
% Computes the shock values for constrained forecasts necessary to keep
% endogenous variables at their constrained paths
%
% INPUTS:
% - cL             [integer]              scalar, number of controlled periods
% - H              [integer]              scalar, number of forecast periods
% - mcValue        [double]               n_controlled_vars*cL array, paths for constrained variables
% - shocks         [double]               n_controlled_vars*cL array, shock values draws (with zeros for controlled_varexo)
% - forcs          [double]               n_endovars*(H+1) matrix of endogenous variables storing the inital condition
% - T              [double]               n_endovars*n_endovars array, transition matrix of the state equation.
% - R              [double]               n_endovars*n_exo array, matrix relating the endogenous variables to the innovations in the state equation.
% - mv             [logical]              n_controlled_exo*n_endovars array, indicator selecting constrained endogenous variables
% - mu             [logical]              n_controlled_vars*nexo array, indicator selecting controlled exogenous variables
%
% OUTPUTS:
% - forcs          [double]               n_endovars*(H+1) array, forecasted endogenous variables
% - e              [double]               nexo*H array, exogenous variables
%
% ALGORITHM:
%
%   Relies on state-space form:
%
%       yₜ = T  yₜ₋₁ + R  εₜ
%
%   Both yₜ, the vector of endogenous variables, and εₜ are split up into controlled
%   and  uncontrolled ones, and we assume, without loss of generality, that the
%   constrained endogenous variables and the controlled shocks come first :
%
%      ⎧ y₁ₜ ⎫    ⎧ T₁₁  T₁₂ ⎫ ⎧ y₁ₜ₋₁ ⎫   ⎧ R₁₁  R₁₂ ⎫ ⎧ ε₁ₜ ⎫
%      ⎩ y₂ₜ ⎭ =  ⎩ T₂₁  T₂₂ ⎭ ⎩ y₂ₜ₋₁ ⎭ + ⎩ R₂₁  R₂₂ ⎭ ⎩ ε₂ₜ ⎭
%
%   where matrices T and R are partitioned consistently with the
%   vectors of endogenous variables and innovations. Provided that matrix
%   R₁₁ is square and full rank (a necessary condition is that the
%   number of free endogenous variables matches the number of free innovations),
%   given y₁ₜ, ε₂ₜ and yₜ₋₁ the first block of equations can be solved for ε₁ₜ:
%
%      ε₁ₜ = R₁₁ \ ( y₁ₜ - T₁₁y₁ₜ₋₁ - T₁₂y₂ₜ₋₁  - R₁₂ε₂ₜ )
%
%   and y₂ₜ can be updated by evaluating the second block of equations:
%
%      y₂ₜ = T₂₁y₁ₜ₋₁ + T₂₂y₂ₜ₋₁ +  R₂₁ε₁ₜ + R₂₂ε₂ₜ
%
%   By iterating over these two blocks of equations, we can build a forecast for
%   all the endogenous variables in the system conditional on paths for a subset of the
%   endogenous variables. This exercise is replicated by drawing different
%   sequences of free innovations. The result is a predictive distribution for
%   the uncontrolled endogenous variables, y₂ₜ, that Dynare will use to report
%   confidence bands around the point conditional forecast.
%   is used for forecasting

% Copyright © 2006-2022 Dynare Team
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

    if cL
        e = zeros(size(mcValue,1),cL);
        for t = 1:cL % Loop over the two blocks of equations
            k = find(isfinite(mcValue(:,t))); % missing conditional values are indicated by NaN
            e(k,t) = inv(mv(k,:)*R*mu(:,k))*(mcValue(k,t)-mv(k,:)*T*forcs(:,t)-mv(k,:)*R*shocks(:,t));
            forcs(:,t+1) = T*forcs(:,t)+R*(mu(:,k)*e(k,t)+shocks(:,t));
        end
    end
    for t = cL+1:H
        forcs(:,t+1) = T*forcs(:,t)+R*shocks(:,t);
    end