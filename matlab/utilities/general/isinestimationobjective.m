function b = isinestimationobjective()

% Returns true iff the caller function is called by a routine defining the objective of an estimation algorithm.

% Copyright © 2014-2020 Dynare Team
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

b = false;

tmp = struct2cell(dbstack);
list_of_previously_called_routines = transpose(tmp(2,:));

if ~isempty(strmatch('dsge_likelihood', list_of_previously_called_routines, 'exact')) || ...
        ~isempty(strmatch('dsge_var_likelihood', list_of_previously_called_routines, 'exact')) || ...
        ~isempty(strmatch('non_linear_dsge_likelihood', list_of_previously_called_routines, 'exact'))
    b = true;
end
