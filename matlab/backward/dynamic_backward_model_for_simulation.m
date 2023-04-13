function [r, J] = dynamic_backward_model_for_simulation(z, dynamic_resid, dynamic_g1, ylag, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr)

% Dynamic routine's wrapper used by dynare_solve for simulating backward models

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

endo_nbr = length(z);

% Build y vector to be passed to the dynamic model.
y = [ ylag; z; NaN(endo_nbr, 1) ];

[r, T_order, T] = dynamic_resid(y, x, params, steady_state);

if nargout>1
    Jacobian = dynamic_g1(y, x, params, steady_state, sparse_rowval, ...
                          sparse_colval, sparse_colptr, T_order, T);
    J = Jacobian(:, endo_nbr+(1:endo_nbr));
end
