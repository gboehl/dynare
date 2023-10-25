function print_info(info, noprint, options_)
% print_info(info, noprint, options_)
% Prints error messages
%
% INPUTS
%   info              [double]     vector returned by resol.m
%   noprint           [integer]    equal to 0 if the error message has to be printed.
%   options_          [structure]  Dynare options
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2005-2023 Dynare Team
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
if ~noprint
    message = get_error_message(info, options_);
    error(message);
end