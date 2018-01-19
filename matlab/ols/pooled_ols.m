function varargout = pooled_ols(ds, param_common, param_regex, overlapping_dates, eqtags)
% function pooled_ols(ds, param_common, param_regex, overlapping_dates, eqtags)
% Run Pooled OLS
% Apply parameter values found to corresponding parameter values in the
% other blocks of the model
%
% INPUTS
%   ds                  [dseries]  data to use in estimation
%   param_common        [cellstr]  List of values to insert into param_regex,
%                                  e.g. country codes {'FR', 'DE', 'IT'}
%   param_regex         [cellstr]  Where '*' should be replaced by the first
%                                  value in param_common
%   overlapping_dates   [bool]     if the dates across the equations should
%                                  overlap
%   eqtags              [cellstr]  names of equation tags to estimate. If empty,
%                                  estimate all equations
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   dynare must be run with the option: json=compute

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

global M_ oo_

%% Check input arguments
assert(~isempty(ds) && isdseries(ds), 'The first argument must be a dseries');

if isempty(param_common) && isempty(param_regex)
    disp('Performing OLS instead of Pooled OLS...')
    if nargin < 6
        dyn_ols(ds);
    else
        dyn_ols(ds, {}, eqtags);
    end
    return;
end
assert(~isempty(param_common) && iscellstr(param_common), 'The second argument must be a cellstr');
assert(~isempty(param_regex) && iscellstr(param_regex), 'The third argument must be a cellstr');

if nargin < 4
    overlapping_dates = false;
else
    assert(islogical(overlapping_dates) && length(overlapping_dates) == 1, 'The fourth argument must be a bool');
end

%% Read JSON
jsonfile = [M_.fname '_original.json'];
if exist(jsonfile, 'file') ~= 2
    error('Could not find %s! Please use the json=compute option (See the Dynare invocation section in the reference manual).', jsonfile);
end

jsonmodel = loadjson(jsonfile);
jsonmodel = jsonmodel.model;
if nargin < 5
    [lhs, rhs, lineno] = getEquationsByTags(jsonmodel);
else
    [lhs, rhs, lineno] = getEquationsByTags(jsonmodel, 'name', eqtags);
end

%% Replace parameter names in equations
country_name = param_common{1};
regexcountries = ['(' strjoin(param_common(2:end),'|') ')'];
param_regex_idx = false(length(param_regex), 1);
for i = 1:length(param_regex)
    splitp = strsplit(param_regex{i}, '*');
    assert(length(splitp) >= 2);
    rhstmp = regexprep(rhs, ...
        strjoin(splitp, regexcountries), ...
        strjoin(splitp, country_name));
    if length(intersect(rhs, rhstmp)) ~= length(rhs)
        rhs = rhstmp;
        param_regex_idx(i) = true;
    end
end
param_regex = param_regex(param_regex_idx);

st = dbstack(1);
save_structure_name = 'pooled_ols';
if strcmp(st(1).name, 'pooled_fgls')
    varargout{1} = param_regex;
    save_structure_name = 'pooled_fgls';
end

%% Find parameters and variable names in every equation & Setup estimation matrices
[X, Y, startdates, enddates, startidxs, residnames, pbeta, vars, surpidxs, surconstrainedparams] = ...
    pooled_sur_common(ds, lhs, rhs, lineno);

if overlapping_dates
    maxfp = max([startdates{:}]);
    minlp = min([enddates{:}]);
    nobs = minlp - maxfp;
    newY = zeros(nobs*length(lhs), 1);
    newX = zeros(nobs*length(lhs), columns(X));
    newstartidxs = zeros(size(startidxs));
    newstartidxs(1) = 1;
    for i = 1:length(lhs)
        if i == length(lhs)
            yds = dseries(Y(startidxs(i):end), startdates{i});
            xds = dseries(X(startidxs(i):end, :), startdates{i});
        else
            yds = dseries(Y(startidxs(i):startidxs(i+1)-1), startdates{i});
            xds = dseries(X(startidxs(i):startidxs(i+1)-1, :), startdates{i});
        end
        newY(newstartidxs(i):newstartidxs(i) + nobs, 1) = yds(maxfp:minlp).data;
        newX(newstartidxs(i):newstartidxs(i) + nobs, :) = xds(maxfp:minlp, :).data;
        if i ~= length(lhs)
            newstartidxs(i+1) = newstartidxs(i) + nobs + 1;
        end
    end
    Y = newY;
    X = newX;
    startidxs = newstartidxs;
    oo_.(save_structure_name).sample_range = maxfp:minlp;
    oo_.(save_structure_name).residnames = residnames;
    oo_.(save_structure_name).Y = Y;
    oo_.(save_structure_name).X = X;
    oo_.(save_structure_name).pbeta = pbeta;
    oo_.(save_structure_name).country_name = country_name;
end

%% Estimation
% Estimated Parameters
[q, r] = qr(X, 0);
oo_.(save_structure_name).beta = r\(q'*Y);

if strcmp(st(1).name, 'pooled_fgls')
    return
end

% Assign parameter values back to parameters using param_regex & param_common
regexcountries = ['(' strjoin(param_common(1:end),'|') ')'];
assigned_idxs = false(size(pbeta));
for i = 1:length(param_regex)
    beta_idx = strcmp(pbeta, strrep(param_regex{i}, '*', country_name));
    assigned_idxs = assigned_idxs | beta_idx;
    value = oo_.(save_structure_name).beta(beta_idx);
    assert(~isempty(value));
    M_.params(~cellfun(@isempty, regexp(M_.param_names, ...
        strrep(param_regex{i}, '*', regexcountries)))) = value;
end
idxs = find(assigned_idxs == 0);
values = oo_.(save_structure_name).beta(idxs);
names = pbeta(idxs);
assert(length(values) == length(names));
for i = 1:length(idxs)
    M_.params(strcmp(M_.param_names, names{i})) = values(i);
end

residuals = Y - X * oo_.(save_structure_name).beta;
for i = 1:length(lhs)
    if i == length(lhs)
        oo_.(save_structure_name).resid.(residnames{i}{:}) = residuals(startidxs(i):end);
    else
        oo_.(save_structure_name).resid.(residnames{i}{:}) = residuals(startidxs(i):startidxs(i+1)-1);
    end
    oo_.(save_structure_name).varcovar.(['eq' num2str(i)]) = oo_.(save_structure_name).resid.(residnames{i}{:})*oo_.(save_structure_name).resid.(residnames{i}{:})';
    idx = find(strcmp(residnames{i}{:}, M_.exo_names));
    M_.Sigma_e(idx, idx) = var(oo_.(save_structure_name).resid.(residnames{i}{:}));
end
end
