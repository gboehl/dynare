function info = isquare(A)

% Returns true iff A is a square matrix.
%
% INPUTS
% - A       [double]     matrix.
%
% OUTPUTS
% - info    [logical]

% Copyright © 2013-2018 Dynare Team
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

info = false;
if ismatrix(A) && isequal(size(A, 1), size(A, 2))
    info = true;
end