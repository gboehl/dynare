function g4_unfolded = unfold_g4(g4, ny)
% Given the 4th order derivatives stored in a sparse matrix and without
% symmetric elements (as returned by the static/dynamic files) and the number
% of (static or dynamic) variables in the jacobian, returns
% an unfolded version of the same matrix (i.e. with symmetric elements).

% Copyright (C) 2019 Dynare Team
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

[i, j, v] = find(g4);

i_unfolded = [];
j_unfolded = [];
v_unfolded = [];

for k = 1:length(v)
    l1 = rem(j(k)-1, ny);
    j2 = floor((j(k)-1)/ny);
    l2 = rem(j2, ny);
    j3 = floor((j(k)-1)/ny^2);
    l3 = rem(j3,ny);
    l4 = floor(j3/ny);

    p = unique(perms([l1 l2 l3 l4]), 'rows');
    np = rows(p);

    i_unfolded = [i_unfolded; repmat(i(k), np, 1)];
    j_unfolded = [j_unfolded; 1 + p(:,1) + ny*(p(:,2) + ny*(p(:,3) + ny*p(:,4)))];
    v_unfolded = [v_unfolded; repmat(v(k), np, 1)];
end

g4_unfolded = sparse(i_unfolded, j_unfolded, v_unfolded, size(g4, 1), size(g4, 2));
