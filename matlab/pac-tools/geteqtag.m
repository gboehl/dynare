function eqtag = geteqtag(eqname, pacname, Model)

% Returns the equation tag, as defined by the preprocessor,
% associated to a pac equation.
%
% INPUTS
% - eqname     [char]    1×n array, name of the PAC equation.
% - pacname    [char]    1×m array, name of the PAC model.
%
% OUTPUTS
% - eqtag      [char]    equation tag associated to eqname.

% Copyright © 2019 Dynare Team
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

eqtag = Model.pac.(pacname).tag_map{strcmp(eqname, Model.pac.(pacname).tag_map(:,1)),2};