function pooled_ols(ds, param_common, param_regex, overlapping_dates, save_structure_name)
% function pooled_ols(ds, param_common, param_regex, overlapping_dates, save_structure_name)
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
%   save_structure_name [string]   Name of structure in oo_ to save results in
%                                  (pooled_ols by default)
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

% Check input arguments
assert(~isempty(ds) && isdseries(ds), 'The first argument must be a dseries');

if isempty(param_common) && isempty(param_regex)
    disp('Performing OLS instead of Pooled OLS...')
    dyn_ols(ds);
    return;
end
assert(~isempty(param_common) && iscellstr(param_common), 'The second argument must be a cellstr');
assert(~isempty(param_regex) && iscellstr(param_regex), 'The third argument must be a cellstr');

if nargin < 4
    overlapping_dates = false;
else
    assert(islogical(overlapping_dates) && length(overlapping_dates) == 1, 'The fourth argument must be a bool');
end

if nargin < 5
    save_structure_name = 'pooled_ols';
else
    assert(ischar(save_structure_name), 'The fifth argument must be a string');
end

%% Read JSON
jsonfile = [M_.fname '_original.json'];
if exist(jsonfile, 'file') ~= 2
    error('Could not find %s! Please use the json=parse option (See the Dynare invocation section in the reference manual).', jsonfile);
end

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
M_exo_names_trim = cellstr(M_.exo_names);
M_endo_exo_names_trim = [cellstr(M_.endo_names); M_exo_names_trim];
M_param_names_trim = cellstr(M_.param_names);
regex = strjoin(M_endo_exo_names_trim(:,1), '|');
mathops = '[\+\*\^\-\/]';
params = cell(length(rhs),1);
vars = cell(length(rhs),1);
pbeta = {};
Y = [];
X = [];
startidxs = zeros(length(lhs), 1);
startdates = cell(length(lhs), 1);
enddates = cell(length(lhs), 1);
residnames = cell(length(lhs), 1);
for i = 1:length(lhs)
    rhs_ = strsplit(rhs{i}, {'+','-','*','/','^','log(','ln(','log10(','exp(','(',')','diff('});
    rhs_(cellfun(@(x) all(isstrprop(x, 'digit')), rhs_)) = [];
    vnames = setdiff(rhs_, M_param_names_trim);
    if ~isempty(regexp(rhs{i}, ...
            ['(' strjoin(vnames, '\\(\\d+\\)|') '\\(\\d+\\))'], ...
            'once'))
        error(['pooled_ols: you cannot have leads in equation on line ' ...
            lineno{i} ': ' lhs{i} ' = ' rhs{i}]);
    end

    % Find parameters and associated variables
    pnames = intersect(rhs_, M_param_names_trim);
    pidxs = zeros(length(pnames), 1);
    vnames = cell(1, length(pnames));
    splitstrings = cell(length(pnames), 1);
    xjdata = dseries;
    for j = 1:length(pnames)
        createdvar = false;
        idx = find(strcmp(pbeta, pnames{j}));
        if isempty(idx)
            pbeta = [pbeta; pnames{j}];
            pidxs(j) = length(pbeta);
        else
            pidxs(j) = idx;
        end

        pregex = [...
            mathops pnames{j} mathops ...
            '|^' pnames{j} mathops ...
            '|' mathops pnames{j} '$' ...
            ];
        [startidx, endidx] = regexp(rhs{i}, pregex, 'start', 'end');
        assert(length(startidx) == 1);
        if rhs{i}(startidx) == '*' && rhs{i}(endidx) == '*'
            vnames{j} = [getStrMoveLeft(rhs{i}(1:startidx-1)) '*' ...
                getStrMoveRight(rhs{i}(endidx+1:end))];
        elseif rhs{i}(startidx) == '*'
            vnames{j} = getStrMoveLeft(rhs{i}(1:startidx-1));
            splitstrings{j} = [vnames{j} '*' pnames{j}];
        elseif rhs{i}(endidx) == '*'
            vnames{j} = getStrMoveRight(rhs{i}(endidx+1:end));
            splitstrings{j} = [pnames{j} '*' vnames{j}];
            if rhs{i}(startidx) == '-'
                vnames{j} = ['-' vnames{j}];
                splitstrings{j} = ['-' splitstrings{j}];
            end
        elseif rhs{i}(startidx) == '+' ...
                || rhs{i}(startidx) == '-' ...
                || rhs{i}(endidx) == '+' ...
                || rhs{i}(endidx) == '-'
            % intercept
            createdvar = true;
            if any(strcmp(M_endo_exo_names_trim, 'intercept'))
                [~, vnames{j}] = fileparts(tempname);
                vnames{j} = ['intercept_' vnames{j}];
                assert(~any(strcmp(M_endo_exo_names_trim, vnames{j})));
            else
                vnames{j} = 'intercept';
            end
            splitstrings{j} = vnames{j};
        else
            error('pooled_ols: Shouldn''t arrive here');
        end
        if createdvar
            xjdatatmp = dseries(ones(ds.nobs, 1), ds.firstdate, vnames{j});
        else
            xjdatatmp = eval(regexprep(vnames{j}, regex, 'ds.$&'));
            xjdatatmp.rename_(vnames{j});
        end
        xjdatatmp.rename_(num2str(j));
        xjdata = [xjdata xjdatatmp];
    end

    lhssub = getRhsToSubFromLhs(ds, rhs{i}, regex, [splitstrings; pnames]);
    
    residuals = intersect(rhs_, cellstr(M_.exo_names));
    justvnames = regexprep(vnames, '\(-\d\)|log|exp|log10|[\(\)]', '');
    justvnames = regexp(justvnames, '[-+]', 'split');
    justvnames = [justvnames{:}];
    for j = 1:length(residuals)
        if any(strcmp(residuals{j}, justvnames))
            residuals{j} = [];
        end
    end
    idx = ~cellfun(@isempty, residuals);
    assert(~isempty(idx), ['No residuals in equation ' num2str(i)]);
    assert(sum(idx) == 1, ['More than one residual in equation ' num2str(i)]);
    residnames{i} = residuals{idx};

    params{i} = pnames;
    vars{i} = [vnames{:}];

    ydata = eval(regexprep(lhs{i}, regex, 'ds.$&'));
    for j = 1:lhssub.vobs
        ydata = ydata - lhssub{j};
    end

    if isempty(xjdata)
        % AR(1) case
        fp = ydata.firstobservedperiod;
        lp = ydata.lastobservedperiod;
        startidxs(i) = length(Y) + 1;
        startdates{i} = fp;
        enddates{i} = lp;
        Y(startidxs(i):startidxs(i)+lp-fp, 1) = ydata(fp:lp).data;
        X(startidxs(i):startidxs(i)+lp-fp, :) = zeros(ydata(fp:lp).nobs, columns(X));
    else
        fp = max(ydata.firstobservedperiod, xjdata.firstobservedperiod);
        lp = min(ydata.lastobservedperiod, xjdata.lastobservedperiod);
        
        startidxs(i) = length(Y) + 1;
        startdates{i} = fp;
        enddates{i} = lp;
        Y(startidxs(i):startidxs(i)+lp-fp, 1) = ydata(fp:lp).data;
        X(startidxs(i):startidxs(i)+lp-fp, pidxs) = xjdata(fp:lp).data;
    end
end

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

st = dbstack(1);
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
    M_.params(~cellfun(@isempty, regexp(M_param_names_trim, ...
        strrep(param_regex{i}, '*', regexcountries)))) = value;
end
idxs = find(assigned_idxs == 0);
values = oo_.(save_structure_name).beta(idxs);
names = pbeta(idxs);
assert(length(values) == length(names));
for i = 1:length(idxs)
    M_.params(strcmp(M_param_names_trim, names{i})) = values(i);
end

residuals = Y - X * oo_.(save_structure_name).beta;
for i = 1:length(lhs)
    if i == length(lhs)
        oo_.(save_structure_name).resid.(residnames{i}) = residuals(startidxs(i):end);
    else
        oo_.(save_structure_name).resid.(residnames{i}) = residuals(startidxs(i):startidxs(i+1)-1);
    end
    oo_.(save_structure_name).varcovar.(['eq' num2str(i)]) = oo_.(save_structure_name).resid.(residnames{i})*oo_.(save_structure_name).resid.(residnames{i})';
    idx = find(strcmp(residnames{i}, M_exo_names_trim));
    M_.Sigma_e(idx, idx) = var(oo_.(save_structure_name).resid.(residnames{i}));
end
end
