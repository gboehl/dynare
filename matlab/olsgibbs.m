function ds = olsgibbs(ds, eqtag, BetaPriorExpectation, BetaPriorVariance, s2, nu, ndraws, discarddraws, thin, fitted_names_dict)
%function ds = olsgibbs(ds, eqtag, BetaPriorExpectation, BetaPriorVariance, s2, nu, ndraws, discarddraws, thin, fitted_names_dict)
% Implements Gibbs Sampling for univariate linear model.
%
% INPUTS
% - ds                          [dseries]    dataset.
% - eqtag                       [string]     name of equation tag to estimate.
% - BetaPriorExpectation        [double]     vector with n elements, prior expectation of β.
% - BetaPriorVariance           [double]     n*n matrix, prior variance of β.
% - s2                          [double]     scalar, first hyperparameter for h.
% - nu                          [integer]    scalar, second hyperparameter for h.
% - ndraws                      [integer]    scalar, total number of draws (Gibbs sampling)
% - discarddraws                [integer]    scalar, number of draws to be discarded.
% - thin                        [integer]    scalar, if thin == N, save every Nth draw (default is 1).
% - fitted_names_dict           [cell]       Nx2 or Nx3 cell array to be used in naming fitted
%                                            values; first column is the equation tag,
%                                            second column is the name of the
%                                            associated fitted value, third column
%                                            (if it exists) is the function name of
%                                            the transformation to perform on the
%                                            fitted value.
%
% OUTPUTS
% - ds                          [dseries]    dataset updated with fitted value
%
% SPECIAL REQUIREMENTS
%   dynare must have been run with the option: json=compute

% Copyright (C) 2018-2019 Dynare Team
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

global M_ oo_ options_

%% Check input
if nargin < 7 || nargin > 10
    error('Incorrect number of arguments passed to olsgibbs')
end

if isempty(ds) || ~isdseries(ds)
    error('The 1st argument must be a dseries')
end

if ~ischar(eqtag)
    error('The 2nd argument must be a string')
end

if ~isvector(BetaPriorExpectation)
    error('The 3rd argument must be a vector')
else
    if ~isempty(BetaPriorExpectation)
        BetaPriorExpectation = transpose(BetaPriorExpectation(:));
    end
end

if ~ismatrix(BetaPriorVariance) || length(BetaPriorExpectation)~=length(BetaPriorVariance)
    error('The 4th argument (BetaPriorVariance) must be a square matrix with the same dimension as the third argument (BetaPriorExpectation)')
else
    warning('off', 'MATLAB:singularMatrix')
    BetaInversePriorVariance = eye(length(BetaPriorVariance))/BetaPriorVariance;
    warning('on', 'MATLAB:singularMatrix')
end

if ~isreal(s2)
    error('The 5th argument (s2) must be a double')
end

if ~isint(nu)
    error('The 6th argument (nu) must be an integer')
end

if ~isint(ndraws)
    error('The 7th argument (ndraws) must be an integer')
end

if nargin <= 7
    discarddraws = 0;
else
    if ~isint(discarddraws)
        error('The 8th argument (discardeddraws), if provided, must be an integer')
    else
        if discarddraws >= ndraws
            error('The 8th argument (discardeddraws) must be smaller than the 7th argument (ndraws)')
        end
    end
end

if nargin <= 8
    thin = 1;
else
    if ~isint(thin)
        error('The 9th argument, must be an integer')
    end
end

if nargin <= 9
    fitted_names_dict = {};
else
    if ~isempty(fitted_names_dict) && ...
            (~iscell(fitted_names_dict) || ...
            (size(fitted_names_dict, 2) < 2 || size(fitted_names_dict, 2) > 3))
        error('The 10th argument must be an Nx2 or Nx3 cell array');
    end
end

%% Parse equation
[Y, lhssub, X, fp, lp] = common_parsing(ds, get_ast({eqtag}), true);
lhsname = Y{1}.name;
Y = Y{1}.data;
X = X{1};
fp = fp{1};
lp = lp{1};
pnames = X.name;
N = size(X.data, 1);
X = X.data;

%% Estimation (see Koop, Gary. Bayesian Econometrics. 2003. Chapter 4.2)
PosteriorDegreesOfFreedom = N + nu;
n = length(pnames);
assert(n==length(BetaPriorExpectation), 'the length prior mean for beta must be the same as the number of parameters in the equation to be estimated.');
h = 1.0/s2*nu; % Initialize h to the prior expectation.

periods = 1;
linee = 1;

% Posterior Simulation
oo_.olsgibbs.(eqtag).draws = zeros(floor((ndraws-discarddraws)/thin), n+3);
for i=1:discarddraws
    % Set conditional distribution of β
    InverseConditionalPoseriorVariance = BetaInversePriorVariance + h*(X'*X);
    cholConditionalPosteriorVariance = chol(InverseConditionalPoseriorVariance\eye(n), 'upper');
    ConditionalPosteriorExpectation = (BetaPriorExpectation*BetaInversePriorVariance + h*(Y'*X))/InverseConditionalPoseriorVariance;
    % Draw beta | Y, h
    beta = rand_multivariate_normal(ConditionalPosteriorExpectation, cholConditionalPosteriorVariance, n);
    % draw h | Y, beta
    resids = Y - X*transpose(beta);
    s2_ = (resids'*resids + nu*s2)/PosteriorDegreesOfFreedom;
    h = gamrnd(PosteriorDegreesOfFreedom/2.0, 2.0/(PosteriorDegreesOfFreedom*s2_));
end

for i = discarddraws+1:ndraws
    % Set conditional distribution of β
    InverseConditionalPoseriorVariance = BetaInversePriorVariance + h*(X'*X);
    cholConditionalPosteriorVariance = chol(InverseConditionalPoseriorVariance\eye(n), 'upper');
    ConditionalPosteriorExpectation = (BetaPriorExpectation*BetaInversePriorVariance + h*(Y'*X))/InverseConditionalPoseriorVariance;
    % Draw beta | Y, h
    beta = rand_multivariate_normal(ConditionalPosteriorExpectation, cholConditionalPosteriorVariance, n);
    % draw h | Y, beta
    resids = Y - X*transpose(beta);
    s2_ = (resids'*resids + nu*s2)/PosteriorDegreesOfFreedom;
    h = gamrnd(PosteriorDegreesOfFreedom/2.0, 2.0/(PosteriorDegreesOfFreedom*s2_));
    R2 = 1-var(resids)/var(Y);
    if isequal(periods, thin)
        oo_.olsgibbs.(eqtag).draws(linee, 1:n) = beta;
        oo_.olsgibbs.(eqtag).draws(linee, n+1) = h;
        oo_.olsgibbs.(eqtag).draws(linee, n+2) = s2_;
        oo_.olsgibbs.(eqtag).draws(linee, n+3) = R2;
        periods = 1;
        linee = linee+1;
    else
        periods = periods+1;
    end
end

%% Save posterior moments.
oo_.olsgibbs.(eqtag).posterior.mean.beta = mean(oo_.olsgibbs.(eqtag).draws(:,1:n))';
oo_.olsgibbs.(eqtag).posterior.mean.h = mean(oo_.olsgibbs.(eqtag).draws(:,n+1));
oo_.olsgibbs.(eqtag).posterior.variance.beta = cov(oo_.olsgibbs.(eqtag).draws(:,1:n));
oo_.olsgibbs.(eqtag).posterior.variance.h = var(oo_.olsgibbs.(eqtag).draws(:,n+1));
oo_.olsgibbs.(eqtag).s2 = mean(oo_.olsgibbs.(eqtag).draws(:,n+2));
oo_.olsgibbs.(eqtag).R2 = mean(oo_.olsgibbs.(eqtag).draws(:,n+3));

% Yhat
idx = 0;
yhatname = [eqtag '_olsgibbs_FIT'];
if ~isempty(fitted_names_dict)
    idx = strcmp(fitted_names_dict(:,1), eqtag);
    if any(idx)
        yhatname = fitted_names_dict{idx, 2};
    end
end
oo_.olsgibbs.(eqtag).Yhat = dseries(X*oo_.olsgibbs.(eqtag).posterior.mean.beta, fp, yhatname) + lhssub{1};

% Apply correcting function for Yhat if it was passed
if any(idx) ...
        && length(fitted_names_dict(idx, :)) == 3 ...
        && ~isempty(fitted_names_dict{idx, 3})
    oo_.olsgibbs.(eqtag).Yhat = ...
        feval(fitted_names_dict{idx, 3}, oo_.olsgibbs.(eqtag).Yhat);
end
ds = [ds oo_.olsgibbs.(eqtag).Yhat];

% Compute and save posterior densities.
for i=1:n
    xx = oo_.olsgibbs.(eqtag).draws(:,i);
    nn = length(xx);
    bandwidth = mh_optimal_bandwidth(xx, nn, 0, 'gaussian');
    [x, f] = kernel_density_estimate(xx, 512, nn, bandwidth,'gaussian');
    oo_.olsgibbs.(eqtag).posterior.density.(pnames{i}) = [x, f];
end

% Update model's parameters with posterior mean.
idxs = zeros(length(pnames), 1);
for j = 1:length(pnames)
    idxs(j) = find(strcmp(M_.param_names, pnames{j}));
    M_.params(idxs(j)) = oo_.olsgibbs.(eqtag).posterior.mean.beta(j);
end

% Write .inc file
write_param_init_inc_file('olsgibbs', eqtag, idxs, oo_.olsgibbs.(eqtag).posterior.mean.beta);

%% Print Output
if ~options_.noprint
    ttitle = ['Bayesian estimation (with Gibbs sampling) of equation ''' eqtag ''''];
    preamble = {['Dependent Variable: ' lhsname{:}], ...
                sprintf('No. Independent Variables: %d', size(X,2)), ...
                sprintf('Observations: %d from %s to %s\n', size(X,1), fp.char, lp.char)};
    afterward = {sprintf('s^2: %f', oo_.olsgibbs.(eqtag).s2), sprintf('R^2: %f', oo_.olsgibbs.(eqtag).R2)};
    dyn_table(ttitle, preamble, afterward, pnames, {'Posterior mean', 'Posterior std.'}, 4, [oo_.olsgibbs.(eqtag).posterior.mean.beta, sqrt(diag(oo_.olsgibbs.(eqtag).posterior.variance.beta))]);
end
end
