function list = list_of_parameters_calibrated_as_NaN(M_)

% The name of the function is explicit enough...
%
% INPUTS
%   M_    [structure]   Description of the (simulated or estimated) model.
%
% OUTPUTS
%   list  [cell]        Each element is the name (row char array) of parameter without value.
%
% ALGORITHM
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2010-2018 Dynare Team
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

idx  = find(isnan(M_.params));
nnn  = length(idx);

if nnn
    list = M_.param_names(idx);
else
    list = {};
end