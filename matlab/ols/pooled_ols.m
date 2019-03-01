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
%   return arguments common to pooled_fgls only if called from pooled_fgls
%   
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

global M_ oo_

%% Check input arguments
if isempty(ds) || ~isdseries(ds)
    error('The first argument must be a dseries');
end

if nargin < 5
    eqtags = {};
end

if isempty(param_common) && isempty(param_regex)
    disp('Performing OLS instead of Pooled OLS...')
    dyn_ols(ds, {}, eqtags);
    return
end

assert(~isempty(param_common) && iscellstr(param_common), 'The second argument must be a cellstr');
assert(~isempty(param_regex) && iscellstr(param_regex), 'The third argument must be a cellstr');

if nargin < 4
    overlapping_dates = false;
else
    assert(islogical(overlapping_dates) && length(overlapping_dates) == 1, 'The fourth argument must be a bool');
end

%% Get Equation(s)
ast = get_ast(eqtags);
neqs = length(ast);

%% Replace parameter names in equations
country_name = param_common{1};
regexcountries = ['(' strjoin(param_common(2:end),'|') ')'];
ast = replace_parameters(ast, country_name, regexcountries, param_regex);

%% Find parameters and variable names in every equation & Setup estimation matrices
[Y, lhssub, X, ~, ~, residnames] = common_parsing(ds, ast, overlapping_dates);
clear ast
nobs = zeros(length(Y), 1);
nobs(1) = Y{1}.nobs;
fp = Y{1}.firstobservedperiod;
for i = 2:length(Y)
    if Y{i}.firstobservedperiod < fp
        fp = Y{i}.firstobservedperiod;
    end
    nobs(i) = Y{i}.nobs;
end
[Y, ~, X] = put_in_sur_form(Y, lhssub, X);

%% Handle FGLS
st = dbstack(1);
if strcmp(st(1).name, 'pooled_fgls')
    save_structure_name = 'pooled_fgls';
    % Pass vars back to pooled_fgls
    varargout{1} = Y.data;
    varargout{2} = X.data;
    varargout{3} = X.name;
    varargout{4} = residnames;
    varargout{5}= country_name;
else
    save_structure_name = 'pooled_ols';
end

%% Estimation
% Estimated Parameters
[q, r] = qr(X.data, 0);
oo_.(save_structure_name).beta = r\(q'*Y.data);

if strcmp(st(1).name, 'pooled_fgls')
    return
end
clear save_structure_name;

% Assign parameter values back to parameters using param_regex & param_common
regexcountries = ['(' strjoin(param_common(1:end),'|') ')'];
assigned_idxs = false(size(X.name));
incidxs = [];
for i = 1:length(param_regex)
    beta_idx = strcmp(X.name, strrep(param_regex{i}, '*', country_name));
    assigned_idxs = assigned_idxs | beta_idx;
    value = oo_.pooled_ols.beta(beta_idx);
    if isempty(eqtags)
        assert(~isempty(value));
    end
    if ~isempty(value)
        idxs = find(~cellfun(@isempty, regexp(M_.param_names, strrep(param_regex{i}, '*', regexcountries))))';
        incidxs = [incidxs idxs];
        M_.params(idxs) = value;
    end
end
idxs = find(assigned_idxs == 0);
values = oo_.pooled_ols.beta(idxs);
names = X.name(idxs);
assert(length(values) == length(names));
for i = 1:length(idxs)
    incidxs = [incidxs find(strcmp(M_.param_names, names{i}))];
    M_.params(incidxs(end)) = values(i);
end

% Write .inc file
write_param_init_inc_file('pooled_ols', M_.fname, incidxs, M_.params(incidxs));

residuals = Y.data - X.data * oo_.pooled_ols.beta;
for i = 1:neqs
    if i == 1
        oo_.pooled_ols.resid.(residnames{i}) = residuals(1:nobs(1));
    elseif i == neqs
        oo_.pooled_ols.resid.(residnames{i}) = residuals(sum(nobs(1:i-1))+1:end);
    else
        oo_.pooled_ols.resid.(residnames{i}) = residuals(sum(nobs(1:i-1))+1:sum(nobs(1:i)));
    end
    oo_.pooled_ols.varcovar.(['eq' num2str(i)]) = oo_.pooled_ols.resid.(residnames{i})*oo_.pooled_ols.resid.(residnames{i})';
    idx = find(strcmp(residnames{i}, M_.exo_names));
    M_.Sigma_e(idx, idx) = var(oo_.pooled_ols.resid.(residnames{i}));
end
end

function ast = replace_parameters(ast, country_name, regexcountries, param_regex)
for i = 1:length(ast)
    ast{i}.AST = replace_parameters_recursive(ast{i}.AST, country_name, regexcountries, param_regex);
end
end

function node = replace_parameters_recursive(node, country_name, regexcountries, param_regex)
if strcmp(node.node_type, 'VariableNode')
    if strcmp(node.type, 'parameter')
        for i = 1:length(param_regex)
            splitp = strsplit(param_regex{i}, '*');
            assert(length(splitp) >= 2);
            tmp = regexprep(node.name, strjoin(splitp, regexcountries), strjoin(splitp, country_name));
            if ~strcmp(tmp, node.name)
                node.name = tmp;
                return
            end
        end
    end
elseif strcmp(node.node_type, 'UnaryOpNode')
    node.arg = replace_parameters_recursive(node.arg, country_name, regexcountries, param_regex);
elseif strcmp(node.node_type, 'BinaryOpNode')
    node.arg1 = replace_parameters_recursive(node.arg1, country_name, regexcountries, param_regex);
    node.arg2 = replace_parameters_recursive(node.arg2, country_name, regexcountries, param_regex);
end
end
