function eqnumber = get_equation_number_by_tag(eqname)

% Translates an equation name into an equation number.
%
% INPUTS
% - eqname   [string]    Name of the equation.
%
% OUTPUTS
% - eqnumber [integer]   Equation number.

% Copyright (C) 2018 Dynare Team
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

eqnumber = strmatch(eqname, M_.equations_tags(strmatch('name', M_.equations_tags(:,2), 'exact'), 3), 'exact');

if isempty(eqnumber), eqnumber = 0; end