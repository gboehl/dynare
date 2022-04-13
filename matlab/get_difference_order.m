function n = get_difference_order(var)

% Returns true iff endogenous variable `var` is a variable in difference.
%
% INPUTS
% - var   [string, integer]  Variable name or index in M_.endo_names.
%
% OUTPUTS
% - boo   [logical]          true/false.

% Copyright © 2018 Dynare Team
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

% Get index of the auxiliary
id = get_aux_variable_id(var);

% Set default difference order
n = 0;

% If var is not an auxiliary variable it cannot be a first difference.
if ~id, return, end

while id
    if M_.aux_vars(id).type==8
        % diff auxiliary variable.
        n = n+1;
        v = M_.aux_vars(id).orig_index;
        if v<=M_.orig_endo_nbr
            id = 0;
        else
            id = get_aux_variable_id(v);
        end
    elseif M_.aux_vars(id).type==9
        % lagged diff auxiliary variable
        v = M_.aux_vars(id).orig_index;
        id = get_aux_variable_id(v);
    else
        break
    end
end