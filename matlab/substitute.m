function v = substitute(v, i, x)

% Substitute a scalar in a vector.
%
% INPUTS
% - v     [double]      m×1 vector
% - i     [integer]     scalar, index for the scalar to be replaced
% - x     [double]      scalar or 1×n vector.
%
% OUTPUTS
% - v     [double]      m×1 vector or m×n matrix (with substituted value(s))
%
% REMARKS
% If x is a vector with n elements, then n substitutions are performed (returning n updated vectors in a matrix with n columns)
%
% EXAMPLES
% >> v = ones(2,1);
% >> substitude(v, 1, 0)
%
% ans =                                                                                                                                                                        %
%
%    0
%
%    1
%
% >> substitute(v, 1, [3 4])
%
% ans =
%
%    3     4
%    1     1

% Copyright © 2023 Dynare Team
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

assert(isvector(v), 'First input argument must be a vector.')
assert(isvector(x), 'Last input argument must be a scalar or a vector.')
assert(isscalar(i) && isint(i) && i>0 && i<=length(v), 'Second input argument must be a scalar integer')

if length(x)==1
    v(i) = x;
else
    v = repmat(v, 1, length(x));
    v(i,:) = x;
end
