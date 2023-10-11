function print_bytecode_dynamic_model()
% function print_bytecode_dynamic_model()
% print the model and jacobian from the bytecode format for the dynamic model
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

global M_ options_ oo_

if options_.bytecode
    z = repmat(oo_.steady_state, 1 , M_.maximum_lead + M_.maximum_lag + 1);
    zx = repmat([oo_.exo_steady_state; oo_.exo_det_steady_state]', M_.maximum_lead + M_.maximum_lag + 1, 1);

    if options_.block
        bytecode('print', 'dynamic', 'block_decomposed', M_, options_, z, zx, M_.params, oo_.steady_state, 1);
    else
        bytecode('print', 'dynamic', M_, options_, z, zx, M_.params, oo_.steady_state, 1);
    end
else
    disp('You have to use bytecode option in model command to use print_bytecode_dynamic_model');
end
