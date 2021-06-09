function A = isconst(B)

% Returns an array of logicals, i-th element is `true` if column i of B is constant, `false` otherwise.
%
% INPUTS:
% - B        [double]    n×m matrix.
%
% OUTPUTS
% - A        [logical]   1×m vector.

% Copyright © 2008-2019 Dynare Team
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
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>

A = false(1, size(B, 2));

for i = 1:size(B, 2)
    if all(abs(B(2:end,i)-B(1:end-1,i))<1e-12)
        A(i) = true;
    end
end