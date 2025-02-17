function x = get_param_by_name(pname)
% function x = get_param_by_name(pname)
% returns the value of a parameter identified by its name
%
% INPUTS:
%   pname:  parameter name
%
% OUTPUTS
%   x:      parameter value
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2006-2017 Dynare Team
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

i = strmatch(pname,M_.param_names,'exact');

if isempty(i)
    error('Can''t find parameter %s', pname)
end

x = M_.params(i);
