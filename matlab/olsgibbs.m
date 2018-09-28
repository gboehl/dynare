function olsgibbs(ds, eqtag, beta0, V, s2, nu, ndraws, discarddraws, thin)
%function olsgibbs(ds, eqtag, beta, V, s2, nu, ndraws, discarddraws, thin)
% Implements Gibbs Samipling for SUR
%
% INPUTS
%   ds           [dseries]    data
%   eqtag        [string]     name of equation tag to estimate.
%   beta0        [vector]     prior mean of beta
%   V            [matrix]     prior covariance of beta
%   s2           [double]     prior mean of h
%   nu           [int]        degrees of freedom of h
%   ndraws       [int]        number of draws
%   discarddraws [int]        number of draws to discard
%   thin         [int]        if thin == N, save every Nth draw
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2018 Dynare Team
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
if nargin < 7 || nargin > 9
    error('Incorrect number of arguments passed to olsgibbs')
end
if ~isdseries(ds)
    error('The 1st argument must be a dseries')
end
if ~ischar(eqtag)
    error('The 2nd argument must be a string')
end
if ~isvector(beta0)
    error('The 3rd argument must be a vector')
else
    if isrow(beta0)
        beta0 = beta0';
    end
end
if ~ismatrix(V) || size(beta0,1) ~= size(V,1) || size(beta0,1) ~= size(V,2)
    error('The 4th argument must be a square matrix with the same dimension as beta0')
else
    V = inv(V);
end
if ~isreal(s2)
    error('The 5th argument must be a double')
end
if ~isint(nu)
    error('The 6th argument must be an integer')
end
if ~isint(ndraws)
    error('The 7th argument must be an integer')
end
if nargin == 7
    discarddraws = 0;
else
    if ~isint(discarddraws)
        error('The 8th argument, if provided, must be an integer')
    end
end
if nargin == 8
    thin = 1;
else
    if ~isint(thin)
        error('The 9th argument, if provided, must be an integer')
    end
end

%% Let dyn_ols parse
[N, pnames, X, Y] = dyn_ols(ds, {}, {eqtag});

%% Estimation
%  Notation from: Koop, Gary. Bayesian Econometrics. 2003. Chapter 4.2

% Setup
nubar = N + nu;
nparams = length(pnames);
assert(nparams == size(beta0,1), ['the length prior mean for beta must '...
    'be the same as the number of parameters in the equation to be estimated.']);

h = rand;
thinidx = 1;
drawidx = 1;

% Posterior Simulation
oo_.olsgibbs.(eqtag).betadraws = zeros(floor((ndraws-discarddraws)/thin), nparams);
for i = 1:ndraws
    % draw beta | Y, h
    Vbar = inv(V + h*(X'*X));
    betabar = Vbar*(V*beta0 + h*X'*Y);
    beta = rand_multivariate_normal(betabar', chol(Vbar), nparams)';

    % draw h | Y, beta
    resid = Y - X*beta;
    s2bar = (resid'*resid + nu*s2)/nubar;
    h = gamma_specification(s2bar, nubar);
    if i > discarddraws
        if thinidx == thin
            oo_.olsgibbs.(eqtag).betadraws(drawidx, :) = beta';
            thinidx = 1;
            drawidx = drawidx + 1;
        else
            thinidx = thinidx + 1;
        end
    end
end

%% Save parameter values
oo_.olsgibbs.(eqtag).beta = (sum(oo_.olsgibbs.(eqtag).betadraws)/rows(oo_.olsgibbs.(eqtag).betadraws))';
for j = 1:length(pnames)
    oo_.olsgibbs.(eqtag).param_idxs(j) = find(strcmp(M_.param_names, pnames{j}));
    M_.params(oo_.ols.(eqtag).param_idxs(j)) = oo_.ols.(eqtag).beta(j);
end
oo_.olsgibbs.(eqtag).s2 = s2bar;
oo_.olsgibbs.(eqtag).h = h;

%% Print Output
if ~options_.noprint
    ttitle = ['Gibbs sampling of linear model in equation ''' eqtag ''''];
    afterward = {sprintf('s^2: %f', oo_.olsgibbs.(eqtag).s2)};
    dyn_table(ttitle, {}, afterward, pnames, ...
        {'Parameter Value'}, 4, oo_.olsgibbs.(eqtag).beta);
    
    % Plot
    figure
    nrows = 5;
    ncols = floor(nparams/nrows);
    if mod(nparams, nrows) ~= 0
        ncols = ncols + 1;
    end
    for j = 1:length(pnames)
        M_.params(strmatch(pnames{j}, M_.param_names, 'exact')) = oo_.olsgibbs.(eqtag).beta(j);
        subplot(nrows, ncols, j)
        histogram(oo_.olsgibbs.(eqtag).betadraws(:, j))
        hc = histcounts(oo_.olsgibbs.(eqtag).betadraws(:, j));
        line([oo_.olsgibbs.(eqtag).beta(j) oo_.olsgibbs.(eqtag).beta(j)], [min(hc) max(hc)], 'Color', 'red');
        title(pnames{j}, 'Interpreter', 'none')
    end
end
end

