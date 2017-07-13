function surgibbs(ds, A, ndraws, varargin)
%function sur(ds)
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
assert(nargin == 1 || nargin == 3, 'Incorrect number of arguments passed to sur');

jsonfile = [M_.fname '_original.json'];
if exist(jsonfile, 'file') ~= 2
    error('Could not find %s! Please use the json option (See the Dynare invocation section in the reference manual).', jsonfile);
end

%% Get Equations
jsonmodel = loadjson(jsonfile);
jsonmodel = jsonmodel.model;
[lhs, rhs, lineno] = getEquationsByTags(jsonmodel, varargin{:});

m = length(lhs);
if m <= 1
    error('SUR estimation requires the selection of at least two equations')
end

%% Construct regression matrices
Y = dseries();
Xi =  cell(m, 1);
pnamesall = [];
vwlagsall = [];
for i = 1:m
    Y = [Y ds{lhs{i}}];

    rhs_ = strsplit(rhs{i}, {'+','-','*','/','^','log(','exp(','(',')'});
    rhs_(cellfun(@(x) all(isstrprop(x, 'digit')), rhs_)) = [];
    vnames = setdiff(rhs_, cellstr(M_.param_names));
    regexprnoleads = cell2mat(strcat('(', vnames, {'\(\d+\))|'}));
    if ~isempty(regexp(rhs{i}, regexprnoleads(1:end-1), 'match'))
        error(['olseqs: you cannot have leads in equation on line ' ...
            lineno{i} ': ' lhs{i} ' = ' rhs{i}]);
    end
    regexpr = cell2mat(strcat('(', vnames, {'\(-\d+\))|'}));
    vwlags = regexp(rhs{i}, regexpr(1:end-1), 'match');

    % Find parameters
    pnames = cell(1, length(vwlags));
    for j = 1:length(vwlags)
        regexmatch = regexp(rhs{i}, ['(\w*\*?)?' strrep(strrep(vwlags{j}, '(', '\('), ')', '\)') '(\*?\w*)?'], 'match');
        regexmatch = strsplit(regexmatch{:}, '*');
        assert(length(regexmatch) == 2);
        if strcmp(vwlags{j}, regexmatch{1})
            pnames{j} = regexmatch{2};
        else
            pnames{j} = regexmatch{1};
        end
    end
    pnamesall = [pnamesall pnames];
    vwlagsall = [vwlagsall vwlags];
    Xi{i} = cellfun(@eval, strcat('ds.', vwlags), 'UniformOutput', false);
end

fp = Y.firstobservedperiod;
lp = Y.lastobservedperiod;
for i = 1:m
    X = dseries();
    for j = 1:length(Xi{i})
        X = [X dseries(Xi{i}{j}.data, Xi{i}{j}.dates, ['V' num2str(i) num2str(j)])];
    end
    Xi{i} = X;
    fp = max(fp, X.firstobservedperiod);
    lp = min(lp, X.lastobservedperiod);
end
Y = Y(fp:lp).data(:);
X = [];
for i = 1:m
    Xi{i} = Xi{i}(fp:lp).data;
    ind = size(X);
    X(ind(1)+1:ind(1)+size(Xi{i}, 1), ind(2)+1:ind(2)+size(Xi{i},2)) = Xi{i};
end

%% Estimation
nobs = length(fp:lp);
nvars = size(X, 2);

[q, r] = qr(X, 0);
resid = Y - X * (r\(q'*Y));
resid = reshape(resid, nobs, m);
S = resid'*resid/nobs;
tmp = kron(inv(S), eye(nobs));
beta0 = (X'*tmp*X)\X'*tmp*Y;
beta = beta0;

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
