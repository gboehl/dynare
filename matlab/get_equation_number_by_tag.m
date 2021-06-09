function eqnumber = get_equation_number_by_tag(eqname, DynareModel)

% Translates an equation name into an equation number.
%
% INPUTS
% - eqname        [char]     1×n array, name of the equation.
% - DynareModel   [struct]   Structure describing the model, aka M_.
%
% OUTPUTS
% - eqnumber      [integer]  Equation number.

% Copyright © 2018-2020 Dynare Team
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

eqnumber = strmatch(eqname, DynareModel.equations_tags(strmatch('name', DynareModel.equations_tags(:,2), 'exact'), 3), 'exact');

if isempty(eqnumber), eqnumber = 0; end