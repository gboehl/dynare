function [Y, lhssub, X, startdates, enddates, residnames] = common_parsing(ds, ast, overlapping_dates, param_names)
%function [Y, lhssub, X, startdates, enddates, residnames] = common_parsing(ds, ast, overlapping_dates, param_names)
%
% Code common to sur.m and pooled_ols.m
%
% INPUTS
%   ds                [dseries]         dataset
%   ast               [cell array]      JSON representation of abstract syntax tree
%   overlapping_dates [boolean]         if true, dates are same across equations
%
% OUTPUTS
%   Y                 [cell array]      dependent variables
%   lhssub            [cell array]      RHS subtracted from LHS
%   X                 [cell array]      regressors
%   startdates        [cell array]      first observed period for each
%                                       equation
%   enddates          [cell array]      last observed period for each
%                                       equation
%   startidxs         [vector]          rows corresponding to each
%                                       equation's observations
%   residnames        [cell array]      name of residual in each equation
%   param_names       [cellstr or cell] list of parameters to estimate (if
%                                       empty, estimate all)
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2019 Dynare Team
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

global M_

%% Initialize variables
neqs = length(ast);
Y = cell(neqs, 1);
lhssub = cell(neqs, 1);
X = cell(neqs, 1);
startdates = cell(neqs, 1);
enddates = cell(neqs, 1);
residnames = cell(neqs, 1);

%% Loop over equations
for i = 1:neqs
    [Y{i}, lhssub{i}, X{i}, residnames{i}, startdates{i}, enddates{i}] = ...
        parse_ols_style_equation(ds, ast{i});
end

if overlapping_dates
    maxfp = max([startdates{:}]);
    minlp = min([enddates{:}]);
    for i = 1:neqs
        Y{i} = Y{i}(maxfp:minlp);
        if ~isempty(X{i})
            X{i} = X{i}(maxfp:minlp);
        end
        if ~isempty(lhssub{i})
            lhssub{i} = lhssub{i}(maxfp:minlp);
        end
        startdates{i} = maxfp;
        enddates{i} = minlp;
    end
end

if ~isempty(param_names)
    if iscell(param_names) && ~iscellstr(param_names)
        assert(length(param_names) == neqs, 'error w param_names arg');
    else
        pn = param_names;
        pn_found_cellstr = false(length(param_names), 1);
    end
    for i = 1:neqs
        if iscell(param_names) && ~iscellstr(param_names)
            pn = param_names{i};
            pn_found_cellstr_i = false(length(pn), 1);
        end
        names = X{i}.name;
        newlhssub = dseries();
        for j = 1:length(names)
            idx = find(strcmp(names{j}, pn));
            if isempty(idx)
                pval = M_.params(strcmp(names{j}, M_.param_names));
                if isnan(pval) || isinf(pval)
                    error(['Could not find param init value for ' names{j}])
                end
                newlhssub = newlhssub + pval * X{i}.(names{j});
                X{i} = X{i}.remove(names{j});
            else
                if iscell(param_names) && ~iscellstr(param_names)
                    pn_found_cellstr_i = true;
                else
                    pn_found_cellstr(idx) = true;
                end
            end
        end
        Y{i} = Y{i} - newlhssub;
        lhssub{i} = lhssub{i} + newlhssub;
        if iscell(param_names) && ~iscellstr(param_names)
            if ~all(pn_found_cellstr_i)
                error(['parameters specified in param_names (' ...
                    strjoin(pn(~pn_found_cellstr_i), ', ') ...
                    ') were not found in eq ' num2str(i) ' to be estimated'])
            end
        end
    end
    if ~(iscell(param_names) && ~iscellstr(param_names))
        if ~all(pn_found_cellstr)
            error(['parameters specified in param_names (' ...
                strjoin(pn(~pn_found_cellstr), ', ') ...
                ') were not found in the equation(s) to be estimated'])
        end
    end
end
end
