function nds = evaluate(ds, eqtag)

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

global M_

% Get equation
[LHS, RHS] = get_lhs_and_rhs(eqtag, M_, true);

% Parse equation and return list of parameters, endogenous and exogenous variables.
[pnames, enames, xnames] = get_variables_and_parameters_in_equation(LHS, RHS, M_);

% Load parameter values.
commands = sprintf('%s = %s;', pnames{1}, num2str(M_.params(strcmp(pnames{1},M_.param_names)), 16));
for i=2:length(pnames)
    commands = sprintf('%s %s = %s;', commands, pnames{i}, ...
                       num2str(M_.params(strcmp(pnames{i},M_.param_names)), 16));
end
eval(commands)


% Substitute endogenous variable x with ds.x
enames = unique(enames);
for i=1:length(enames)
    if ismember(enames{i}, ds.name)
        RHS = regexprep(RHS, sprintf('\\<(%s)\\>', enames{i}), sprintf('ds.%s', enames{i}));
    else
        error('Endogenous variable %s is unknown in dseries objet.', enames{i})
    end
end

% Substitute exogenous variable x with ds.x, except if
if ~isfield(M_, 'simulation_exo_names')
    M_.simulation_exo_names = M_.exo_names;
end
xnames = unique(xnames);
for i=1:length(xnames)
    if ismember(xnames{i}, M_.simulation_exo_names)
        if ismember(xnames{i}, ds.name)
            RHS = regexprep(RHS, sprintf('\\<(%s)\\>', xnames{i}), sprintf('ds.%s', xnames{i}));
        else
            RHS = regexprep(RHS, sprintf('\\<(%s)\\>', xnames{i}), '0');
            warning('Exogenous variable %s is unknown in dseries objet. Assign zero value.', xnames{i})
        end
    else
        RHS = regexprep(RHS, sprintf('(\\ *)(+)(\\ *)%s', xnames{i}), '');
    end
end

nds = eval(RHS);