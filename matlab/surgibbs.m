function surgibbs(ds, param_names, beta0, A, ndraws, discarddraws, thin, eqtags)
%function surgibbs(ds, param_names, beta0, A, ndraws, discarddraws, thin, eqtags)
% Implements Gibbs Samipling for SUR
%
% INPUTS
%   ds           [dseries]    data
%   param_names  [cellstr]    list of parameters to estimate
%   beta0        [vector]     prior values (in order of param_names)
%   A            [matrix]     prior distribution variance (in order of
%                             param_names)
%   ndraws       [int]        number of draws
%   discarddraws [int]        number of draws to discard
%   thin         [int]        if thin == N, save every Nth draw
%   eqtags       [cellstr]    names of equation tags to estimate. If empty,
%                             estimate all equations
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   dynare must have been run with the option: json=compute

% Copyright (C) 2017-2019 Dynare Team
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
% Bayesian Analysis Volume 5, Number 1, pp. 65-96.

global M_ oo_ options_

%% Check input
assert(nargin >= 5 && nargin <= 8, 'Incorrect number of arguments passed to surgibbs');
assert(isdseries(ds), 'The 1st argument must be a dseries');
assert(iscellstr(param_names), 'The 2nd argument must be a cellstr');
assert(isvector(beta0) && length(beta0) == length(param_names), ...
    'The 3rd argument must be a vector with the same length as param_names and the same ');
if isrow(beta0)
    beta0 = beta0';
end
assert(ismatrix(A) && all(all((A == A'))) && length(beta0) == size(A, 2), ...
    'The 4th argument must be a symmetric matrix with the same dimension as beta0');
assert(isint(ndraws), 'The 5th argument must be an integer');
if nargin == 5
    discarddraws = 0;
else
    assert(isint(discarddraws), 'The 6th argument, if provided, must be an integer');
end
if nargin == 6
    thin = 1;
else
    assert(isint(thin), 'The 7th argument, if provided, must be an integer');
end

%% Estimation
%  Notation from:
%  Ando, Tomohiro and Zellner, Arnold. Hierarchical Bayesian Analysis of
%  the Seemingly Unrelated Regression and Simultaneous Equations Models
%  Using a Combination of Direct Monte Carlo and Importance Sampling
%  Techniques. Bayesian Analysis. 2010. pp 67-70.
if nargin == 8
    [nobs, X, Y, m] = sur(ds, param_names, eqtags);
else
    [nobs, X, Y, m] = sur(ds, param_names);
end

beta = beta0;
A = inv(A);
thinidx = 1;
drawidx = 1;
nparams = length(param_names);
oo_.surgibbs.betadraws = zeros(floor((ndraws-discarddraws)/thin), nparams);
if ~options_.noprint
    disp('surgibbs: estimating, please wait...')
end
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
    beta = rand_multivariate_normal(betabar', chol(Omegabar), nparams)';
    if i > discarddraws
        if thinidx == thin
            oo_.surgibbs.betadraws(drawidx, :) = beta';
            thinidx = 1;
            drawidx = drawidx + 1;
        else
            thinidx = thinidx + 1;
        end
    end
end

% save parameter values
oo_.surgibbs.beta = (sum(oo_.surgibbs.betadraws)/rows(oo_.surgibbs.betadraws))';

incidxs = zeros(length(param_names), 1);
for i = 1:length(param_names)
    incidxs(i) = strmatch(param_names{i}, M_.param_names, 'exact');
    M_.params(incidxs(i)) = oo_.surgibbs.beta(i);
end

% Write .inc file
write_param_init_inc_file('surgibbs', M_.fname, incidxs, oo_.surgibbs.beta);

%% Print Output
if ~options_.noprint
    dyn_table('Gibbs Sampling on SUR', {}, {}, param_names, ...
        {'Parameter Value'}, 4, oo_.surgibbs.beta);
end

%% Plot
if ~options_.nograph
    figure
    nrows = 5;
    ncols = floor(nparams/nrows);
    if mod(nparams, nrows) ~= 0
        ncols = ncols + 1;
    end
    for j = 1:length(param_names)
        subplot(nrows, ncols, j)
        histogram(oo_.surgibbs.betadraws(:, j))
        hc = histcounts(oo_.surgibbs.betadraws(:, j));
        line([oo_.surgibbs.beta(j) oo_.surgibbs.beta(j)], [min(hc) max(hc)], 'Color', 'red');
        title(param_names{j}, 'Interpreter', 'none')
    end
end
end

