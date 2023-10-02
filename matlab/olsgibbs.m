function ds = olsgibbs(ds, eqtag, BetaPriorExpectation, BetaPriorVariance, s2, nu, ndraws, discarddraws, thin, fitted_names_dict, model_name, param_names, ds_range)
%function ds = olsgibbs(ds, eqtag, BetaPriorExpectation, BetaPriorVariance, s2, nu, ndraws, discarddraws, thin, fitted_names_dict, model_name, param_names, ds_range)
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
% - model_name                  [string]     name to use in oo_ and inc file
% - param_names                 [cellstr]    list of parameters to estimate (if
%                                            empty, estimate all)
% - ds_range                    [dates]      range of dates to use in estimation
%
% OUTPUTS
% - ds                          [dseries]    dataset updated with fitted value
%
% SPECIAL REQUIREMENTS
%   dynare must have been run with the option: json=compute

% Copyright © 2018-2023 Dynare Team
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

global M_ oo_ options_

%% Check input
if nargin < 7 || nargin > 13
    error('Incorrect number of arguments')
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

if nargin <= 10
    model_name = eqtag;
else
    if ~isvarname(model_name)
        error('The 11th argument must be a valid string');
    end
end

if nargin <= 11
    param_names = {};
else
    if ~isempty(param_names) && ~iscellstr(param_names)
        error('The 12th argument, if provided, must be a cellstr')
    end
end

if nargin <= 12
    ds_range = ds.dates;
else
    if isempty(ds_range)
        ds_range = ds.dates;
    else
        if ds_range(1) < ds.firstdate || ds_range(end) > lastdate(ds)
            error('There is a problem with the 13th argument: the date range does not correspond to that of the dseries')
        end
    end
end

%% Parse equation
[Y, lhssub, X, fp, lp] = common_parsing(ds(ds_range), get_ast({eqtag}), true, param_names);
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
oo_.olsgibbs.(model_name).draws = zeros(floor((ndraws-discarddraws)/thin), n+3);
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

hh_fig = dyn_waitbar(0,'Please wait. Gibbs sampler...');
set(hh_fig,'Name','Olsgibbs estimation.');
for i = discarddraws+1:ndraws
    if ~mod(i,100)
        dyn_waitbar((i-discarddraws)/(ndraws-discarddraws),hh_fig,'Please wait. Gibbs sampler...');
    end
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
        oo_.olsgibbs.(model_name).draws(linee, 1:n) = beta;
        oo_.olsgibbs.(model_name).draws(linee, n+1) = h;
        oo_.olsgibbs.(model_name).draws(linee, n+2) = s2_;
        oo_.olsgibbs.(model_name).draws(linee, n+3) = R2;
        periods = 1;
        linee = linee+1;
    else
        periods = periods+1;
    end
end
dyn_waitbar_close(hh_fig);

%% Save posterior moments.
oo_.olsgibbs.(model_name).posterior.mean.beta = mean(oo_.olsgibbs.(model_name).draws(:,1:n))';
oo_.olsgibbs.(model_name).posterior.mean.h = mean(oo_.olsgibbs.(model_name).draws(:,n+1));
oo_.olsgibbs.(model_name).posterior.variance.beta = cov(oo_.olsgibbs.(model_name).draws(:,1:n));
oo_.olsgibbs.(model_name).posterior.variance.h = var(oo_.olsgibbs.(model_name).draws(:,n+1));
oo_.olsgibbs.(model_name).s2 = mean(oo_.olsgibbs.(model_name).draws(:,n+2));
oo_.olsgibbs.(model_name).R2 = mean(oo_.olsgibbs.(model_name).draws(:,n+3));

% Yhat
idx = 0;
yhatname = [eqtag '_olsgibbs_FIT'];
if ~isempty(fitted_names_dict)
    idx = strcmp(fitted_names_dict(:,1), eqtag);
    if any(idx)
        yhatname = fitted_names_dict{idx, 2};
    end
end
oo_.olsgibbs.(model_name).Yhat = dseries(X*oo_.olsgibbs.(model_name).posterior.mean.beta, fp, yhatname);
oo_.olsgibbs.(model_name).YhatOrig = oo_.olsgibbs.(model_name).Yhat;
oo_.olsgibbs.(model_name).Yobs = dseries(Y, fp, lhsname);

% Residuals
oo_.olsgibbs.(model_name).resid = Y - oo_.olsgibbs.(model_name).Yhat;

% Apply correcting function for Yhat if it was passed
oo_.olsgibbs.(model_name).Yhat = oo_.olsgibbs.(model_name).Yhat + lhssub{1};
if any(idx) ...
        && length(fitted_names_dict(idx, :)) == 3 ...
        && ~isempty(fitted_names_dict{idx, 3})
    oo_.olsgibbs.(model_name).Yhat = ...
        feval(fitted_names_dict{idx, 3}, oo_.olsgibbs.(model_name).Yhat);
end
ds.(oo_.olsgibbs.(model_name).Yhat.name{:}) = oo_.olsgibbs.(model_name).Yhat;

% Compute and save posterior densities.
for i=1:n
    xx = oo_.olsgibbs.(model_name).draws(:,i);
    nn = length(xx);
    bandwidth = mh_optimal_bandwidth(xx, nn, 0, 'gaussian');
    [x, f] = kernel_density_estimate(xx, 512, nn, bandwidth,'gaussian');
    oo_.olsgibbs.(model_name).posterior.density.(pnames{i}) = [x, f];
end

% Update model's parameters with posterior mean.
idxs = zeros(length(pnames), 1);
for j = 1:length(pnames)
    idxs(j) = find(strcmp(M_.param_names, pnames{j}));
    M_.params(idxs(j)) = oo_.olsgibbs.(model_name).posterior.mean.beta(j);
end
oo_.olsgibbs.(model_name).pnames = pnames;

% Write .inc file
write_param_init_inc_file('olsgibbs', model_name, idxs, oo_.olsgibbs.(model_name).posterior.mean.beta);

%% Print Output
if ~options_.noprint
    ttitle = ['Bayesian estimation (with Gibbs sampling) of equation ''' eqtag ''''];
    preamble = {['Dependent Variable: ' lhsname{:}], ...
                sprintf('No. Independent Variables: %d', size(X,2)), ...
                sprintf('Observations: %d from %s to %s\n', size(X,1), fp.char, lp.char)};
    afterward = {sprintf('s^2: %f', oo_.olsgibbs.(model_name).s2), sprintf('R^2: %f', oo_.olsgibbs.(model_name).R2)};
    dyn_table(ttitle, preamble, afterward, pnames, {'Posterior mean', 'Posterior std.'}, 4, [oo_.olsgibbs.(model_name).posterior.mean.beta, sqrt(diag(oo_.olsgibbs.(model_name).posterior.variance.beta))]);
end
end
