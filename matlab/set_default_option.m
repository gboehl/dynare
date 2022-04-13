function options=set_default_option(options,field,default)

% function options=set_default_option(options,field,default)
% Sets the option value
%
% INPUTS
%    options
%    field:   option name
%    default: assigns a value
%
% OUTPUTS
%    options
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2003-2018 Dynare Team
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

if ~isfield(options,field)
    options.(field) = default;
    return
end

if isempty(options.(field))
    options.(field) = default;
    return
end

if ~iscell(options.(field)) && ~isdates(options.(field)) && ~isstruct(options.(field))
    if isnan(options.(field))
        options.(field) = default;
        return
    end
end

% 06/07/03 MJ added ; to eval expression
