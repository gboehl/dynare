function H = build_two_dim_hessian(sparse_indices, g2_v, neq, nvar)
% Creates a 2D Hessian (equations in rows, pairs of variables in columns),
% given the output from the sparse {static,dynamic}_g2.m
%
% – sparse_indices is typically equal to M_.{static,dynamic}_g2_sparse_indices
% – g2_v is the vector of non zero values returned by {static,dynamic}_g2.m
% – neq is the number of equations (equal to number of rows of the output matrix)
% – nvar is the number of variables (the output matrix will have nvar² columns)

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

nnz = size(sparse_indices, 1);

%% The g2_* arrays may be expanded if there are symmetric elements added
g2_i = int32(zeros(nnz, 1));
g2_j = int32(zeros(nnz, 1));

next_sym_idx = nnz + 1; % Index of next symmetric element to be added

for k = 1:length(g2_v)
    eq = sparse_indices(k, 1);
    var1 = sparse_indices(k, 2)-1;
    var2 = sparse_indices(k, 3)-1;

    g2_i(k) = eq;
    g2_j(k) = var1 * nvar + var2 + 1;

    %% Add symmetric elements, which are not included by sparse {static,dynamic}_g2.m
    if var1 ~= var2
        g2_i(next_sym_idx) = eq;
        g2_j(next_sym_idx) = var2 * nvar + var1 + 1;
        g2_v(next_sym_idx) = g2_v(k);
        next_sym_idx = next_sym_idx + 1;
    end
end

%% On MATLAB < R2020a, sparse() does not accept int32 indices
if ~isoctave && matlab_ver_less_than('9.8')
    g2_i = double(g2_i);
    g2_j = double(g2_j);
end

H = sparse(g2_i, g2_j, g2_v, neq, nvar*nvar);
