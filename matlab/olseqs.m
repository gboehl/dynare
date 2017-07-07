function olseqs(ds, varargin)
%function olseqs(ds, varargin)
% Run OLS on chosen model equations
%
% INPUTS
%   ds      [dseries]    data
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
assert(nargin == 1 || nargin == 3, 'Incorrect number of arguments passed to olseqs');

jsonfile = [M_.fname '_original.json'];
if exist(jsonfile, 'file') ~= 2
    error('Could not find %s! Please use the json option (See the Dynare invocation section in the reference manual).', jsonfile);
end

%% Get Equation(s)
jsonmodel = loadjson(jsonfile);
jsonmodel = jsonmodel.model;
[lhs, rhs, lineno] = getEquationsByTags(jsonmodel, varargin{:});

for i = 1:length(lhs)
    %% Construct regression matrices
    Y = ds{lhs{i}}.data;

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

    X = cell2mat(cellfun(@eval, strcat('ds.', vwlags, '.data'), 'UniformOutput', false));

    % Remove all rows that have a NaN
    [row, ~] = find(isnan(X), 1, 'last');
    Y = Y(row+1:end, :);
    X = X(row+1:end, :);
    
    %% Estimation
    % From LeSage, James P. "Applied Econometrics using MATLAB"
    if nargin == 3
        if iscell(varargin{2})
            tagv = varargin{2}{i};
        else
            tagv = varargin{2};
        end
    else
        tagv = ['eqlineno' num2str(lineno{i})];
    end
    [nobs, nvars] = size(X);
    oo_.ols.(tagv).dof = nobs - nvars;

    % Estimated Parameters
    [q, r] = qr(X, 0);
    xpxi = (r'*r)\eye(nvars);
    oo_.ols.(tagv).beta = r\(q'*Y);
    for j = 1:length(pnames)
        M_.params(strmatch(pnames{j}, M_.param_names, 'exact')) = oo_.ols.(tagv).beta(j);
    end
    
    % Yhat
    oo_.ols.(tagv).Yhat = X*oo_.ols.(tagv).beta;

    % Residuals
    oo_.ols.(tagv).resid = Y - oo_.ols.(tagv).Yhat;

    %% Calculate statistics
    % Estimate for sigma^2
    SS_res = oo_.ols.(tagv).resid'*oo_.ols.(tagv).resid;
    oo_.ols.(tagv).s2 = SS_res/oo_.ols.(tagv).dof;

    % R^2
    ym = Y - mean(Y);
    SS_tot = ym'*ym;
    oo_.ols.(tagv).R2 = 1 - SS_res/SS_tot;

    % Adjusted R^2
    oo_.ols.(tagv).adjR2 = oo_.ols.(tagv).R2 - (1 - oo_.ols.(tagv).R2)*nvars/(oo_.ols.(tagv).dof-1);

    % Durbin-Watson
    ediff = oo_.ols.(tagv).resid(2:nobs) - oo_.ols.(tagv).resid(1:nobs-1);
    oo_.ols.(tagv).dw = (ediff'*ediff)/SS_res;

    % Standard Error
    oo_.ols.(tagv).stderr = sqrt(oo_.ols.(tagv).s2*diag(xpxi));

    % T-Stat
    oo_.ols.(tagv).tstat = oo_.ols.(tagv).beta./oo_.ols.(tagv).stderr;
    
    %% Print Output
    fprintf('OLS Estimation of equation on line %d', lineno{i});
    if nargin == 3
        fprintf(' [%s = %s]', varargin{1}, tagv);
    end
    fprintf('\nDependent Variable: %s\n', lhs{i});
    fprintf('No. Independent Variables: %d\n', nvars);
    fprintf('Observations: %d\n', nobs);
    maxstrlen = 0;
    for j=1:length(vwlags)
        slen = length(vwlags{j});
        if slen > maxstrlen
            maxstrlen = slen;
        end
    end
    titlespacing = repmat(' ', 1, 4 + maxstrlen + 4) ;
    fprintf('%sCoefficients    t-statistic      Std. Error\n', titlespacing);
    fprintf('%s____________    ____________    ____________\n\n', titlespacing);
    format = ['    %-' num2str(maxstrlen) 's'];
    for j = 1:length(vwlags)
        fprintf(format, vwlags{j});
        fprintf('%12.5f    %12.5f    %12.5f\n', ...
            oo_.ols.(tagv).beta(j), ...
            oo_.ols.(tagv).tstat(j), ...
            oo_.ols.(tagv).stderr(j));
    end
    fprintf('\nR^2: %f\n', oo_.ols.(tagv).R2);
    fprintf('R^2 Adjusted: %f\n', oo_.ols.(tagv).adjR2);
    fprintf('s^2: %f\n', oo_.ols.(tagv).s2);
    fprintf('Durbin-Watson: %f\n', oo_.ols.(tagv).dw);
    fprintf('%s\n', repmat('-', 1, 4 + maxstrlen + 4 + 44));
end
end
