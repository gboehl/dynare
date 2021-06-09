function [ds, json] = evaluate(ds, eqtags, firstperiod, lastperiod, json)

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
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

global M_

debug = false;

if ischar(eqtags)
    eqtags = {eqtags};
end

list_of_expression_tokens = {'+', '-', '*', '/', '^', ...
        'exp(', 'log(', 'sqrt(', 'abs(', 'sign(', ...
        'sin(', 'cos(', 'tan(', 'asin(', 'acos(', 'atan(', ...
        'min(', 'max(', ...
        'normcdf(', 'normpdf(', 'erf(', ...
        'diff(', 'adl(', ')'};

if ismember(nargin, [4, 5])
    if isempty(firstperiod) && isempty(lastperiod)
        range = ds.dates(1):ds.dates(end);
    elseif isempty(firstperiod) && ~isempty(lastperiod)
        range = ds.dates(1):lastperiod;
    elseif ~isempty(firstperiod) && isempty(lastperiod)
        range = firstperiod:ds.dates(end);
    else
        range = firstperiod:lastperiod;
    end
elseif isequal(nargin, 3)
    if isempty(firstperiod)
        range = ds.dates(1):ds.dates(end);
    else
        range = firstperiod:ds.dates(end);
    end
elseif isequal(nargin, 2)
    range = ds.dates(1):ds.dates(end);
else
    error('This routine admits 2, 3, 4, or 5 input arguments.')
end

for i=1:length(eqtags)
    % Get equation
    if isequal(i, 1)
        if nargin<5
            [LHS, RHS, json] = get_lhs_and_rhs(eqtags{i}, M_, true);
        else
            [LHS, RHS] = get_lhs_and_rhs(eqtags{i}, M_, true, json);
        end
    else
        [LHS, RHS] = get_lhs_and_rhs(eqtags{i}, M_, true, json);
    end
    % Parse equation and return list of parameters, endogenous and exogenous variables.
    [pnames, enames, xnames] = get_variables_and_parameters_in_equation(LHS, RHS, M_);
    % Load parameter values.
    if ~isempty(pnames)
        commands = sprintf('%s = %s;', pnames{1}, num2str(M_.params(strcmp(pnames{1}, M_.param_names)), 16));
        for j=2:length(pnames)
            commands = sprintf('%s %s = %s;', commands, pnames{j}, num2str(M_.params(strcmp(pnames{j}, M_.param_names)), 16));
        end
        eval(commands)
    end
    % Remove repetitions in enames
    enames = unique(enames);
    % Test if LHS is an endogenous variable
    is_lhs_expression = ~ismember(LHS, enames);
    if is_lhs_expression
        variable = strsplit(LHS, list_of_expression_tokens);
        variable(cellfun(@(x) all(isempty(x)), variable)) = [];
        if length(variable)>1
            error('It is not possible to have an expression with more than one variable on the LHS (%s).', LHS)
        else
            if isequal(LHS, sprintf('log(%s)', variable{1}))
                transform = {'exp'};
            elseif isequal(LHS, sprintf('diff(%s)', variable{1}))
                transform = {'cumsum'};
            elseif isequal(LHS, sprintf('diff(log(%s))', variable{1}))
                transform = {'cumsum', 'exp'};
            elseif isequal(LHS, sprintf('diff(diff(%s))', variable{1}))
                transform = {'cumsum', 'cumsum'};
            elseif isequal(LHS, sprintf('diff(diff(log(%s)))', variable{1}))
                transform = {'cumsum', 'cumsum', 'exp'};
            else
                error('Cannot proceed with provided LHS (%s in %s)', LHS, eqtags{i})
            end
            lhs = variable{1};
        end
    else
        lhs = LHS;
        transform = {};
    end
    % Throw an error if the equation is dynamic.
    if exactcontains(RHS, lhs)
        error('RHS cannot contain LHS variable (%s in %s)', lhs, eqtags{i})
    end
    % Substitute endogenous variable x with ds.
    for j=1:length(enames)
        if ismember(enames{j}, ds.name)
            RHS = exactstrrep(RHS, enames{j}, sprintf('ds(range).%s', enames{j}));
        else
            RHS = exactstrrep(RHS, sprintf('(%s\\((\\-)*\\d\\)|%s)', enames{j}, enames{j}), '0');
            if debug
                warning off backtrace
                warning('Endogenous variable %s is unknown in dseries objet. Assign zero value.', enames{j})
                warning on backtrace
            end
        end
    end
    % Substitute exogenous variable x with ds.x, except if
    if ~isfield(M_, 'simulation_exo_names')
        M_.simulation_exo_names = M_.exo_names;
    end
    xnames = unique(xnames);
    for j=1:length(xnames)
        if ismember(xnames{j}, M_.simulation_exo_names)
            if ismember(xnames{j}, ds.name)
                RHS = exactstrrep(RHS, xnames{j}, sprintf('ds(range).%s', xnames{j}));
            else
                RHS = exactstrrep(RHS, xnames{j}, '0');
                if debug
                    warning off backtrace
                    warning('Exogenous variable %s is unknown in dseries objet. Assign zero value.', xnames{j})
                    warning on backtrace
                end
            end
        else
            RHS = regexprep(RHS, sprintf('(\\ *)(+)(\\ *)%s', xnames{j}), '');
        end
    end
    if isempty(transform)
        ds{LHS} = eval(RHS);
    else
        tmp = eval(RHS);
        switch length(transform)
          case 1
            if isequal(transform{1}, 'cumsum')
                ds{lhs} = cumsum(tmp)+ds{lhs}(range(1)-1).data;
            else
                ds{lhs} = feval(transform{1}, tmp);
            end
          case 2
            if isequal(transform{2}, 'cumsum')
                % Squared first difference.
                t2 = zeros(length(range), 1);
                for t = 1:length(range)
                    t2(t) = 2*ds{lhs}(range(t)-1).data-ds{lhs}(range(t)-2).data+tmp(range(t)).data;
                end
                ds{lhs} = dseries(t2, range(1));
            else
                t2 = zeros(length(range), 1);
                for t = 1:length(range)
                    t1 = feval(transform{2}, log(ds{lhs}(range(t)-1))+tmp(range(t)).data );
                    t2(t) = t1.data;
                end
                ds{lhs} = dseries(t2, range(1));
% $$$                 % The commented version below is more efficient but the discrepancy with what is returned by simulating
% $$$                 % the model is much bigger (see pac/trend-component-28/example4.mod).
% $$$                 tmp = cumsum(tmp)+log(ds{lhs}(range(1)-1).data);
% $$$                 ds{lhs} = feval(transform{2}, tmp);
            end
          case 3
                t2 = zeros(length(range), 1);
                for t = 1:length(range)
                    t2(t) = feval(transform{3}, 2*log(ds{lhs}(range(t)-1).data)-log(ds{lhs}(range(t)-2).data)+tmp(range(t)).data);
                end
                ds{lhs} = dseries(t2, range(1));
          otherwise
            error('More than 3 unary ops. in LHS not implemented.')
        end
    end
end