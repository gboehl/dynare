function ida = get_aux_variable_id(var)

% Returns the index of an auxiliary variable in M_.aux_vars
%
% INPUTS
% - var   [string, integer]  Variable name or index in M_.endo_names.
%
% OUTPUTS
% - ida   [integer]          Corresponding index in M_.aux_vars.

% Copyright (C) 2018-2019 Dynare Team
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

if isempty(var)
    ida = 0;
    return
end

if ischar(var)
    id = find(strcmp(var, M_.endo_names));
    if isempty(id)
        ida = 0;
        return
    else
        var = id;
    end
else
    if ~isint(var) || var>M_.endo_nbr || var<1
        error('Input must be the name of an endogenous variable or an integer between 1 and %s!', num2str(M_.endo_nbr))
    end
end

if var<=M_.orig_endo_nbr
    % var is not an auxiliary variable.
    ida = 0;
else
    ida = find([M_.aux_vars.endo_index]==var);
end