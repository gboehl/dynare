function surgibbs(ds, A, ndraws, discarddraws)
%function surgibbs(ds, ndraws, discarddraws)
% Implements Gibbs Samipling for SUR
%
% INPUTS
%   ds           [dseries]    data
%   A            [matrix]     prior distribution variance
%   ndraws       [int]        number of draws
%   discarddraws [int]        number of draws to discard
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2017-2018 Dynare Team
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

%% The notation that follows comes from Section 2.2 of
% Ando, Tomohiro and Zellner, Arnold. 2010. Hierarchical Bayesian Analysis of the
% Seemingly Unrelated Regression and Simultaneous Equations Models Using a
% Combination of Direct Monte Carlo and Importance Sampling Techniques.
% Bayesian Analysis Volume 5, Number 1, pp. 65?96.

global M_ oo_

%% Check input
assert(nargin == 3 || nargin == 4, 'Incorrect number of arguments passed to surgibbs');
if nargin == 3
    discarddraws = 0;
end

%% Estimation
beta0 = M_.params;
[nobs, nvars, pidxs, beta, X, Y, m] = sur(ds);
pnamesall = cellstr(M_.param_names(pidxs, :));
if any(isnan(beta0))
    beta0 = beta;
else
    beta = beta0;
end
A = inv(A);
oo_.surgibbs.betadraws = zeros(ndraws-discarddraws, nvars);
for i = 1:ndraws
    % Draw Omega, given X, Y, Beta
    resid = reshape(Y - X*beta, nobs, m);
    Omega = rand_inverse_wishart(m, nobs, chol(inv(resid'*resid/nobs)));

    % Draw beta, given X, Y, Omega
    tmp = kron(inv(Omega), eye(nobs));
    tmp1 = X'*tmp*X;
    Omegabar = inv(tmp1 + A);
    betahat = tmp1\X'*tmp*Y;
    betabar = Omegabar*(tmp1*betahat+A*beta0);
    beta = rand_multivariate_normal(betabar', chol(Omegabar), nvars)';
    if i > discarddraws
        oo_.surgibbs.betadraws(i-discarddraws, :) = beta';
    end
end

% save parameter values
oo_.surgibbs.beta = (sum(oo_.surgibbs.betadraws)/(ndraws-discarddraws))';
M_.params(pidxs, 1) = oo_.surgibbs.beta;

%% Print Output
dyn_table('Gibbs Sampling on SUR', {}, {}, pnamesall, ...
    {'Parameter Value'}, 4, oo_.surgibbs.beta);

%% Plot
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
