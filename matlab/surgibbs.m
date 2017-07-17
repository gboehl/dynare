function surgibbs(ds, A, ndraws, varargin)
%function surgibbs(ds, A, ndraws, varargin)
% Implements Gibbs Samipling for SUR
%
% INPUTS
%   ds      [dseries]    data
%   A       [matrix]     prior distribution variance
%   ndraws  [int]        number of draws
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2017 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global M_ oo_

%% Check input
assert(nargin == 3 || nargin == 5, 'Incorrect number of arguments passed to surgibbs');

%% Estimation
[nobs, nvars, pnamesall, beta, X, Y, m] = sur(ds, varargin{:});
beta0 = beta;

oo_.surgibbs.betadraws = zeros(ndraws, nvars);
for i = 1:ndraws
    % Draw S, given X, Y, Beta
    resid = reshape(Y - X*beta, nobs, m);
    Omega = rand_inverse_wishart(m, nobs, (resid'*resid)/nobs);

    % Draw beta, given X, Y, S
    tmp = kron(inv(Omega), eye(nobs));
    tmp1 = X'*tmp*X;
    Omegabar = inv(tmp1 + A);
    betabar = Omegabar*(tmp1*(tmp1\X'*tmp*Y)+A\beta0);
    Sigma_upper_chol = chol(Omegabar, 'upper');
    beta = rand_multivariate_normal(betabar', Sigma_upper_chol, nvars)';
    oo_.surgibbs.betadraws(i, :) = beta';
end

% save parameter values
oo_.surgibbs.beta = (sum(oo_.surgibbs.betadraws)/ndraws)';

% plot
figure
nrows = 5;
ncols = floor(nvars/nrows);
if mod(nvars, nrows) ~= 0
    ncols = ncols + 1;
end
for j = 1:length(pnamesall)
    M_.params(strmatch(pnamesall{j}, M_.param_names, 'exact')) = oo_.surgibbs.beta(j);
    subplot(nrows, ncols, j)
    histogram(oo_.surgibbs.betadraws(:, j))
    hc = histcounts(oo_.surgibbs.betadraws(:, j));
    line([oo_.surgibbs.beta(j) oo_.surgibbs.beta(j)], [min(hc) max(hc)], 'Color', 'red');
    title(pnamesall{j}, 'Interpreter', 'none')
end
