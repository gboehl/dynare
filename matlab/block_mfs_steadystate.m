function [r, g1] = block_mfs_steadystate(y, fh_static, b, y_all, exo, params, T, M)
% Wrapper around the static files, for use with dynare_solve,
% when block option is given to steady.

% Copyright © 2009-2023 Dynare Team
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

y_all(M.block_structure_stat.block(b).variable) = y;

[~,~,r,g1] = fh_static(y_all, exo, params, M.block_structure_stat.block(b).g1_sparse_rowval, ...
                       M.block_structure_stat.block(b).g1_sparse_colval, ...
                       M.block_structure_stat.block(b).g1_sparse_colptr, T);
