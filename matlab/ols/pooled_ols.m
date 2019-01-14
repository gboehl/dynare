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
assert(~isempty(ds) && isdseries(ds), 'The first argument must be a dseries');

if isempty(param_common) && isempty(param_regex)
    disp('Performing OLS instead of Pooled OLS...')
    if nargin < 6
        dyn_ols(ds);
    else
        dyn_ols(ds, {}, eqtags);
    end
    return
end

if nargin < 5
    eqtags = {};
end

assert(~isempty(param_common) && iscellstr(param_common), 'The second argument must be a cellstr');
assert(~isempty(param_regex) && iscellstr(param_regex), 'The third argument must be a cellstr');

if nargin < 4
    overlapping_dates = false;
else
    assert(islogical(overlapping_dates) && length(overlapping_dates) == 1, 'The fourth argument must be a bool');
end


%% Get Equation(s)
[ast, jsonmodel] = get_ast_jsonmodel(eqtags);
neqs = length(jsonmodel);

%% Replace parameter names in equations
country_name = param_common{1};
regexcountries = ['(' strjoin(param_common(2:end),'|') ')'];
ast = replace_parameters(ast, country_name, regexcountries, param_regex);

%% Handle FGLS
st = dbstack(1);
save_structure_name = 'pooled_ols';
if strcmp(st(1).name, 'pooled_fgls')
    varargout{1} = param_regex;
    save_structure_name = 'pooled_fgls';
end

%% Find parameters and variable names in every equation & Setup estimation matrices
[Y, ~, X] = common_parsing(ds, ast, jsonmodel, overlapping_dates);
clear ast jsonmodel;
nobs = Y{1}.nobs;
[Y, X] = put_in_sur_form(Y, X);

%% Save
oo_.(save_structure_name).sample_range = X.firstdate:X.firstdate+nobs;
%oo_.(save_structure_name).residnames = residnames;
oo_.(save_structure_name).Y = Y.data;
oo_.(save_structure_name).X = X.data;
oo_.(save_structure_name).pbeta = X.name;
oo_.(save_structure_name).country_name = country_name;

%% Estimation
% Estimated Parameters
[q, r] = qr(X.data, 0);
oo_.(save_structure_name).beta = r\(q'*Y.data);

if strcmp(st(1).name, 'pooled_fgls')
    return
end

% Assign parameter values back to parameters using param_regex & param_common
regexcountries = ['(' strjoin(param_common(1:end),'|') ')'];
assigned_idxs = false(size(X.name));
for i = 1:length(param_regex)
    beta_idx = strcmp(X.name, strrep(param_regex{i}, '*', country_name));
    assigned_idxs = assigned_idxs | beta_idx;
    value = oo_.(save_structure_name).beta(beta_idx);
    if isempty(eqtags)
        assert(~isempty(value));
    end
    if ~isempty(value)
        M_.params(~cellfun(@isempty, regexp(M_.param_names, ...
            strrep(param_regex{i}, '*', regexcountries)))) = value;
    end
end
idxs = find(assigned_idxs == 0);
values = oo_.(save_structure_name).beta(idxs);
names = X.name(idxs);
assert(length(values) == length(names));
for i = 1:length(idxs)
    M_.params(strcmp(M_.param_names, names{i})) = values(i);
end

residuals = Y.data - X.data * oo_.(save_structure_name).beta;
for i = 1:neqs
    if i == neqs
        oo_.(save_structure_name).resid.(residnames{i}{:}) = residuals(startidxs(i):end);
    else
        oo_.(save_structure_name).resid.(residnames{i}{:}) = residuals(startidxs(i):startidxs(i+1)-1);
    end
    oo_.(save_structure_name).varcovar.(['eq' num2str(i)]) = oo_.(save_structure_name).resid.(residnames{i}{:})*oo_.(save_structure_name).resid.(residnames{i}{:})';
    idx = find(strcmp(residnames{i}{:}, M_.exo_names));
    M_.Sigma_e(idx, idx) = var(oo_.(save_structure_name).resid.(residnames{i}{:}));
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
