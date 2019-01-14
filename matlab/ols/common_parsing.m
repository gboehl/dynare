function [Y, lhssub, X, startdates, enddates, residnames] = common_parsing(ds, ast, jsonmodel, overlapping_dates)
%function [Y, lhssub, X, startdates, enddates, residnames] = common_parsing(ds, ast, jsonmodel, overlapping_dates)
%
% Code common to sur.m and pooled_ols.m
%
% INPUTS
%   ds                   [dseries]     dataset
%   ast                  [cell array]  JSON representation of abstract syntax tree
%   jsonmodel            [cell array]  JSON representation of model block
%   overlapping_dates    [boolean]     if true, dates are same across equations
%
% OUTPUTS
%   Y                    [cell array]  dependent variables
%   lhssub               [cell array]  RHS to subtract from Y
%   X                    [cell array]  regressors
%   startdates           [cell array]  first observed period for each
%                                      equation
%   enddates             [cell array]  last observed period for each
%                                      equation
%   startidxs            [vector]      rows corresponding to each
%                                      equation's observations
%   residnames           [cell array]  name of residual in each equation
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2019 Dynare Team
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


%% Initialize variables
Y = cell(length(ast), 1);
lhssub = cell(length(ast), 1);
X = cell(length(ast), 1);
residnames = cell(length(ast), 1);
startdates = cell(length(ast), 1);
enddates = cell(length(ast), 1);

%% Loop over equations
neqs = length(ast);
for i = 1:neqs
    %% Parse equation i
    [Y{i}, lhssub{i}, X{i}, residnames{i}] = parse_ols_style_equation(ds, ast{i});

    %% Set start and end dates
    [startdates{i}, enddates{i}] = get_ols_start_end_dates(Y{i}, lhssub{i}, X{i}, jsonmodel{i});
    Y{i} = Y{i}(startdates{i}:enddates{i});
    X{i} = X{i}(startdates{i}:enddates{i});
    if ~isempty(lhssub{i})
        lhssub{i} = lhssub{i}(startdates{i}:enddates{i});
    end
end

if overlapping_dates
    maxfp = max([startdates{:}]);
    minlp = min([enddates{:}]);
    for i = 1:neqs
        Y{i} = Y{i}(maxfp:minlp);
        X{i} = X{i}(maxfp:minlp);
        if ~isempty(lhssub{i})
            lhssub{i} = lhssub{i}(maxfp:minlp);
        end
    end
end
end
