function eqnumber = get_equation_number_by_tag(eqname, M_)
% eqnumber = get_equation_number_by_tag(eqname, M_)
% Translates an equation name into an equation number.
%
% INPUTS
% - eqname        [char]     1×n array, name of the equation.
% - M_            [struct]   Structure describing the model
%
% OUTPUTS
% - eqnumber      [integer]  Equation number.

% Copyright © 2018-2023 Dynare Team
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

idx = find(strcmp(M_.equations_tags(:,2), 'name') & ...
           strcmp(M_.equations_tags(:,3), eqname));

if isempty(idx)
    eqnumber = 0;
else
    eqnumber = M_.equations_tags{idx, 1};
end

end
