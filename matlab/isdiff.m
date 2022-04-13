function boo = isdiff(var)

% Returns true iff endogenous variable `var` is a v ariable in difference.
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

ida = get_aux_variable_id(var);
boo = false;

if ~ida, return, end

if ismember(M_.aux_vars(ida).type, [8, 9])
    boo = true;
end