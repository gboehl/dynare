function display_parameter_values

% Copyright © 2024 Dynare Team
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
my_title='Current parameter values:';
labels = M_.param_names;
headers = {'Parameter'; 'Value'};
lh = cellofchararraymaxlength(labels)+2;
options_.noprint=false;
dyntable(options_, my_title, headers, labels, M_.params, lh, 10, 6);
