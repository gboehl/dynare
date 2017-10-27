function pooled_ols(ds, param_common, param_regex)
% function pooled_ols(ds, param_common, param_regex)
% Run Pooled OLS
% Apply parameter values found to corresponding parameter values in the
% other blocks of the model
%
% INPUTS
%   ds            [dseries]      data to use in estimation
%   param_common  [cellstr]      List of values to insert into param_regex,
%                                e.g. country codes {'FR', 'DE', 'IT'}
%   param_regex   [cellstr]      Where '*' should be replaced by the first
%                                value in param_common
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   dynare must be run with the option: json=parse

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

jsonfile = [M_.fname '.json'];
if exist(jsonfile, 'file') ~= 2
    error('Could not find %s! Please use the json=parse option (See the Dynare invocation section in the reference manual).', jsonfile);
end

%% Read JSON
jsonmodel = loadjson(jsonfile);
jsonmodel = jsonmodel.model;
[lhs, rhs, lineno] = getEquationsByTags(jsonmodel);

%% Replace parameter names in equations
country_name = param_common{1};
regexcountries = ['(' strjoin(param_common(2:end),'|') ')'];
for i = 1:length(param_regex)
    splitp = strsplit(param_regex{i}, '*');
    assert(length(splitp) >= 2);
    rhs = regexprep(rhs, ...
        strjoin(splitp, regexcountries), ...
        strjoin(splitp, country_name));
end

%% Find parameters and variable names in every equation & Setup estimation matrices
regexpr1 = ...
    ['(diff\(\w+(\(\W?\w+\))?\))\*$' ...
    '|' '\((\w+(\(\W?\w+\))?(\W?\w+(\(\W?\w+\))?)*)\)\*$' ...
    '|' '(\w+(\(\W?\w+\))?)\*$' ...
    ];

regexpr2 = ...
    ['^\*(diff\(\w+(\(\W?\w+\))?\))' ...
    '|' '^\*\((\w+(\(\W?\w+\))?(\W?\w+(\(\W?\w+\))?)*)' ...
    '|' '^\*(\w+(\(\W?\w+\))?)'
    ];

M_endo_names_trim = cellfun(@strtrim, num2cell(M_.endo_names(:,:),2),  'Uniform', 0);
regex = ['(?<chb>^|(\-|\(|+|\*|\/|\^))(?<var>' ...
    strjoin(M_endo_names_trim, '|') ...
    ')(?<cha>$|(\)\-|\(|+|\*|\/|\^))'];
params = cell(length(rhs),1);
vars = cell(length(rhs),1);
pbeta = {};
Y = [];
X = [];
for i = 1:length(lhs)
    rhs_ = strsplit(rhs{i}, {'+','-','*','/','^','log(','ln(','log10(','exp(','(',')','diff('});
    rhs_(cellfun(@(x) all(isstrprop(x, 'digit')), rhs_)) = [];
    vnames = setdiff(rhs_, cellstr(M_.param_names));
    if ~isempty(regexp(rhs{i}, ...
            ['(' strjoin(vnames, '\\(\\d+\\)|') '\\(\\d+\\))'], ...
            'once'))
        error(['pooled_ols: you cannot have leads in equation on line ' ...
            lineno{i} ': ' lhs{i} ' = ' rhs{i}]);
    end
    assert(numel(intersect(rhs_, cellstr(M_.exo_names))) == 1);

    % Find parameters and associated variables
    pnames = intersect(rhs_, cellstr(M_.param_names));
    pidxs = zeros(length(pnames), 1);
    vnames = cell(1, length(pnames));
    xjdata = dseries;
    for j = 1:length(pnames)
        idx = find(strcmp(pbeta, pnames(j,:)));
        if isempty(idx)
            pbeta = [pbeta; pnames(j,:)];
            pidxs(j) = length(pbeta);
        else
            pidxs(j) = idx;
        end

        rhs_split = strsplit(rhs{i}, pnames{j});
        assert(length(rhs_split) == 2);
        if ~isempty(rhs_split{1}) && rhs_split{1}(end) == '*'
            vnames(j) = regexp(rhs_split{1}, regexpr1, 'tokens');
        elseif ~isempty(rhs_split{2}) && rhs_split{2}(1) == '*'
            vnames(j) = regexp(rhs_split{2}, regexpr2, 'tokens');
        else
            error('pooled_ols: Shouldn''t arrive here');
        end
        xjdatatmp = getdata(ds, regex, vnames{j}{:});
        xjdatatmp.rename_(num2str(j));
        xjdata = [xjdata xjdatatmp];
    end
    params{i} = pnames;
    vars{i} = [vnames{:}];

    ydata = getdata(ds, regex, lhs{i});

    fp = max(ydata.firstobservedperiod, xjdata.firstobservedperiod);
    lp = min(ydata.lastobservedperiod, xjdata.lastobservedperiod);

    Y(length(Y)+1:length(Y)+1+lp-fp, 1) = ydata(fp:lp).data;
    X(size(X,1)+1:size(X,1)+1+lp-fp, pidxs) = xjdata(fp:lp).data;
end

%% Estimation
% Estimated Parameters
[q, r] = qr(X, 0);
oo_.pooled_ols.beta = r\(q'*Y);
oo_.pooled_ols.param_names = pbeta;

% Assign parameter values back to parameters using param_regex & param_common
param_names_trim = cellfun(@strtrim, num2cell(M_.param_names(:,:),2),  'Uniform', 0);
regexcountries = ['(' strjoin(param_common(1:end),'|') ')'];
assigned_idxs = false(size(pbeta));
for i = 1:length(param_regex)
    beta_idx = strcmp(pbeta, strrep(param_regex{i}, '*', country_name));
    assigned_idxs = assigned_idxs | beta_idx;
    value = oo_.pooled_ols.beta(beta_idx);
    assert(~isempty(value));
    M_.params(~cellfun(@isempty, regexp(param_names_trim, ...
        strrep(param_regex{i}, '*', regexcountries)))) = value;
end
idxs = find(assigned_idxs == 0);
values = oo_.pooled_ols.beta(idxs);
names = pbeta(idxs);
assert(length(values) == length(names));
for i = 1:length(idxs)
    M_.params(strcmp(param_names_trim, names{i})) = values(i);
end
end

function retval = getdata(ds, regex, ser)
if strncmp(ser, 'diff', 4)
    ser = ser(6:end-1);
    lagidx = strfind(ser, '(');
    if isempty(lagidx)
        retval = ds{ser} - ds{ser}(-1);
    else
        lag = str2double(ser(lagidx+1:strfind(ser, ')')-1));
        assert(lag < 0);
        ser = ser(1:lagidx-1);
        retval = ds{ser}(lag) - ds{ser}(lag-1);
    end
else
    retval = eval(regexprep(ser, regex, '$<chb>ds.$<var>$<cha>'));
end
end
