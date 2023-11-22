function ds = surgibbs(ds, param_names, beta0, A, ndraws, discarddraws, thin, eqtags, model_name)

% Implements Gibbs Sampling for SUR
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
%   model_name   [string]     name to use in oo_ and inc file
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   dynare must have been run with the option: json=compute
%
% REFERENCES
% - Ando, Tomohiro and Zellner, Arnold. 2010. Hierarchical Bayesian Analysis of the
%   Seemingly Unrelated Regression and Simultaneous Equations Models Using a
%   Combination of Direct Monte Carlo and Importance Sampling Techniques.
%   Bayesian Analysis Volume 5, Number 1, pp. 65-96.

% Copyright © 2017-2023 Dynare Team
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

%
% Check inputs
%

assert(nargin >= 5 && nargin <= 9, 'Incorrect number of arguments passed to surgibbs');
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

if nargin <= 8
    if ~isfield(oo_, 'surgibbs')
        model_name = 'surgibbs_model_number_1';
    else
        model_name = ['surgibbs_model_number_' num2str(length(fieldnames(oo_.surgibbs)) + 1)];
    end
else
    if ~isvarname(model_name)
        error('The 9th argument must be a valid string');
    end
end

%
% Estimation
%

if nargin == 8
    [nobs, X, Y, m, lhssub, fp] = sur(ds, param_names, eqtags);
else
    [nobs, X, Y, m, lhssub, fp] = sur(ds, param_names);
end

oo_.surgibbs.(model_name).dof = nobs;

beta = beta0;
A = inv(A);
thinidx = 1;
drawidx = 1;
nparams = length(param_names);
oo_.surgibbs.(model_name).betadraws = zeros(floor((ndraws-discarddraws)/thin), nparams);
if ~options_.noprint
    disp('surgibbs: estimating, please wait...')
end

hh_fig = dyn_waitbar(0,'Please wait. Gibbs sampler...');
set(hh_fig,'Name','Surgibbs estimation.');
residdraws = zeros(floor((ndraws-discarddraws)/thin), nobs, m);
for i = 1:ndraws
    if ~mod(i,100)
        dyn_waitbar(i/ndraws,hh_fig,'Please wait. Gibbs sampler...');
    end
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
            oo_.surgibbs.(model_name).betadraws(drawidx, 1:nparams) = beta';
            residdraws(drawidx, 1:nobs, 1:m) = resid;
            thinidx = 1;
            drawidx = drawidx + 1;
        else
            thinidx = thinidx + 1;
        end
    end
end
dyn_waitbar_close(hh_fig);

%
% Save results.
%

oo_.surgibbs.(model_name).posterior.mean.beta = (sum(oo_.surgibbs.(model_name).betadraws)/rows(oo_.surgibbs.(model_name).betadraws))';
oo_.surgibbs.(model_name).posterior.variance.beta = cov(oo_.surgibbs.(model_name).betadraws);

% Yhat
oo_.surgibbs.(model_name).Yhat = X*oo_.surgibbs.(model_name).posterior.mean.beta;
oo_.surgibbs.(model_name).YhatOrig = oo_.surgibbs.(model_name).Yhat;
oo_.surgibbs.(model_name).Yobs = Y;

% Residuals
oo_.surgibbs.(model_name).resid = Y - oo_.surgibbs.(model_name).Yhat;

% Correct Yhat reported back to user
oo_.surgibbs.(model_name).Yhat = oo_.surgibbs.(model_name).Yhat + lhssub;
yhatname = [model_name '_FIT'];
ds.(yhatname) = dseries(oo_.surgibbs.(model_name).Yhat,  fp, yhatname);

% Compute and save posterior densities.
for i=1:nparams
    xx = oo_.surgibbs.(model_name).betadraws(:,i);
    nn = length(xx);
    bandwidth = mh_optimal_bandwidth(xx, nn, 0, 'gaussian');
    [x, f] = kernel_density_estimate(xx, 512, nn, bandwidth, 'gaussian');
    oo_.surgibbs.(model_name).posterior.density.(param_names{i}) = [x, f];
end

% Update model1s parameters with posterior mean.
oo_.surgibbs.(model_name).param_idxs = zeros(length(param_names), 1);
for i = 1:length(param_names)
     if ~strcmp(param_names{i}, 'intercept')
         oo_.surgibbs.(model_name).param_idxs(i) = find(strcmp(M_.param_names, param_names{i}));
         M_.params(oo_.surgibbs.(model_name).param_idxs(i)) = oo_.surgibbs.(model_name).posterior.mean.beta(i);
     end
end
oo_.surgibbs.(model_name).pnames = param_names;
oo_.surgibbs.(model_name).neqs = m;

% Estimate for sigma^2
SS_res = oo_.surgibbs.(model_name).resid'*oo_.surgibbs.(model_name).resid;
oo_.surgibbs.(model_name).s2 = SS_res/oo_.surgibbs.(model_name).dof;

% Set appropriate entries in Sigma_e
posterior_mean_resid = reshape((sum(residdraws))/rows(residdraws), nobs, m);
Sigma_e = posterior_mean_resid'*posterior_mean_resid/oo_.surgibbs.(model_name).dof;

% System R² value of McElroy (1977) - formula from Judge et al. (1986, p. 477)
%
% The R² is computed at the posterior mean of the estimated
% parameters. Maybe it would make more sense to compute a posterior
% distribution for this statistic…
oo_.surgibbs.(model_name).R2 = 1 - (oo_.surgibbs.(model_name).resid' * kron(inv(Sigma_e), eye(nobs)) * oo_.surgibbs.(model_name).resid) ...
                                 / (oo_.surgibbs.(model_name).Yobs' * kron(inv(Sigma_e), eye(nobs)-ones(nobs,nobs)/nobs) * oo_.surgibbs.(model_name).Yobs);

% Write .inc file
write_param_init_inc_file('surgibbs', model_name, oo_.surgibbs.(model_name).param_idxs, oo_.surgibbs.(model_name).posterior.mean.beta);

%
% Print Output
%

if ~options_.noprint
    ttitle = 'Gibbs Sampling on SUR';
    preamble = {['Model name: ' model_name], ...
        sprintf('No. Equations: %d', oo_.surgibbs.(model_name).neqs), ...
        sprintf('No. Independent Variables: %d', size(X, 2)), ...
        sprintf('Observations: %d', oo_.surgibbs.(model_name).dof)};

    afterward = {sprintf('s^2: %f', oo_.surgibbs.(model_name).s2), sprintf('R^2: %f', oo_.surgibbs.(model_name).R2)};
    dyn_table(ttitle, preamble, afterward, param_names,...
             {'Posterior mean', 'Posterior std.'}, 4,...
             [oo_.surgibbs.(model_name).posterior.mean.beta, sqrt(diag(oo_.surgibbs.(model_name).posterior.variance.beta))]);
end

%
% Plot
%

% The histogram() function is not implemented in Octave
if ~options_.nograph && ~isoctave
    figure
    nrows = 5;
    ncols = floor(nparams/nrows);
    if mod(nparams, nrows) ~= 0
        ncols = ncols + 1;
    end
    for j = 1:length(param_names)
        subplot(nrows, ncols, j)
        histogram(oo_.surgibbs.(model_name).betadraws(:, j))
        hc = histcounts(oo_.surgibbs.(model_name).betadraws(:, j));
        line([oo_.surgibbs.(model_name).posterior.mean.beta(j) oo_.surgibbs.(model_name).posterior.mean.beta(j)], [min(hc) max(hc)], 'Color', 'red');
        title(param_names{j}, 'Interpreter', 'none')
    end
end
end
