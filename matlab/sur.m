function varargout = sur(ds, varargin)
%function varargout = sur(ds, varargin)
% Run a Seemingly Unrelated Regression on the provided equations
%
% INPUTS
%   ds        [dseries]     data
%
% OUTPUTS
%   varargout [cell array]  contains the common work between sur and
%                           surgibbs
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
xpxi = (r'*r)\eye(nvars);
resid = Y - X * (r\(q'*Y));
resid = reshape(resid, nobs, m);
s2 = resid'*resid/nobs;
tmp = kron(inv(s2), eye(nobs));
beta = (X'*tmp*X)\X'*tmp*Y;

% if called from surgibbs, return common work
st = dbstack(1);
if strcmp(st(1).name, 'surgibbs')
    varargout{1} = nobs;
    varargout{2} = nvars;
    varargout{3} = pnamesall;
    varargout{4} = beta;
    varargout{5} = X;
    varargout{6} = Y;
    varargout{7} = m;
    return
end

oo_.sur.s2 = s2;
oo_.sur.beta = beta;

for j = 1:length(pnamesall)
    M_.params(strmatch(pnamesall{j}, M_.param_names, 'exact')) = oo_.sur.beta(j);
end

% Yhat
oo_.sur.Yhat = X * oo_.sur.beta;

% Residuals
oo_.sur.resid = Y - oo_.sur.Yhat;

%% Calculate statistics
oo_.sur.dof = nobs;

% Estimate for sigma^2
SS_res = oo_.sur.resid'*oo_.sur.resid;
oo_.sur.s2 = SS_res/oo_.sur.dof;

% R^2
ym = Y - mean(Y);
SS_tot = ym'*ym;
oo_.sur.R2 = 1 - SS_res/SS_tot;

% Adjusted R^2
oo_.sur.adjR2 = oo_.sur.R2 - (1 - oo_.sur.R2)*nvars/(oo_.sur.dof-1);

% Durbin-Watson
ediff = oo_.sur.resid(2:nobs) - oo_.sur.resid(1:nobs-1);
oo_.sur.dw = (ediff'*ediff)/SS_res;

% Standard Error
oo_.sur.stderr = sqrt(oo_.sur.s2*diag(xpxi));

% T-Stat
oo_.sur.tstat = oo_.sur.beta./oo_.sur.stderr;

%% Print Output
fprintf('SUR Estimation');
if nargin == 1
    fprintf(' of all equations')
else
    fprintf(' [%s = {', varargin{1});
    for i = 1:length(varargin{2})
        if i ~= 1
            fprintf(', ');
        end
        fprintf('%s', varargin{2}{i});
    end
    fprintf('}]');
end
fprintf('\n    Dependent Variable: %s\n', lhs{i});
fprintf('    No. Independent Variables: %d\n', nvars);
fprintf('    Observations: %d\n', nobs);
maxstrlen = 0;
for j=1:length(vwlagsall)
    slen = length(vwlagsall{j});
    if slen > maxstrlen
        maxstrlen = slen;
    end
end
titlespacing = repmat(' ', 1, 4 + maxstrlen + 4) ;
fprintf('%sCoefficients    t-statistic      Std. Error\n', titlespacing);
fprintf('%s____________    ____________    ____________\n\n', titlespacing);
format = ['    %-' num2str(maxstrlen) 's'];
for j = 1:length(vwlagsall)
    fprintf(format, vwlagsall{j});
    fprintf('%12.5f    %12.5f    %12.5f\n', ...
        oo_.sur.beta(j), ...
        oo_.sur.tstat(j), ...
        oo_.sur.stderr(j));
end
fprintf('\n    R^2: %f\n', oo_.sur.R2);
fprintf('    R^2 Adjusted: %f\n', oo_.sur.adjR2);
fprintf('    s^2: %f\n', oo_.sur.s2);
fprintf('    Durbin-Watson: %f\n', oo_.sur.dw);
fprintf('%s\n\n', repmat('-', 1, 4 + maxstrlen + 4 + 44));
end
