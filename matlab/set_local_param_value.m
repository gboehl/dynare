function M_ = set_local_param_value(pname,value,M_)
% M_ = set_local_param_value(pname,value,M_)

% Copyright © 2017-2023 Dynare Team
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


i = strmatch(pname, M_.param_names, 'exact');

if isempty(i)
    error(['Parameter name ' pname ' doesn''t exist'])
end

M_.params(i) = value;
