function options_=set_option(options_,options_occbin_,fieldname)
% function options_=set_option(options_,options_occbin_,fieldname)
% Set local option for Occbin
%
% Inputs:
% - options_            [structure]     Matlab's structure containing the options
% - options_occbin_     [structure]     Matlab's structure containing Occbin options
% - fieldname           [string]        name of the options field to set 
%
% Outputs:
% - options_            [structure]     Matlab's structure containing the options

% Copyright © 2021 Dynare Team
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

dot_pos = strfind(fieldname,'.');

if isfield(options_occbin_,fieldname(1:dot_pos-1)) && isfield(options_occbin_.(fieldname(1:dot_pos-1)),fieldname(dot_pos+1:end))
    options_.occbin.(fieldname(1:dot_pos-1)).(fieldname(dot_pos+1:end)) = options_occbin_.(fieldname(1:dot_pos-1)).(fieldname(dot_pos+1:end));
end