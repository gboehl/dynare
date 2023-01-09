function print_bytecode_static_model()
% function print_bytecode_static_model()
% print the model and jacobian from the bytecode format for the static model
%
% INPUTS
%   none
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2001-2023 Dynare Team
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
global options_
if options_.bytecode
    if options_.block
        bytecode('print','static','block_decomposed');
    else
        bytecode('print','static');
    end
else
    disp('You have to use bytecode option in model command to use print_bytecode_static_model');
end
