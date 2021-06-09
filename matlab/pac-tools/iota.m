function i = iota(n, idx)

% Returns a selection vector.
%
% INPUTS 
% - n      [integer]   scalar, dimension of the returned vector.
% - idx    [integer]   vector or scalar, non zero entries indices in the returned vector.
%
% OUTPUTS 
% - i      [integer]   n*1 vector. All elements are zero except those specified in idx.

% Copyright (C) 2018 Dynare Team
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

if ~isscalar(n) || ~isvector(idx)
    error('First input has to be a scalar and second input a vector!')
end

if (n-round(n))>eps(0) || any((idx-round(idx))>eps(0))
    error('Inputs must have integer values!')
end

if n<1
    error('First input argument must be a strictly positive integer!')
end

if any(idx>n)
    error('Elements in second argument cannot be greater than the first argument!')
end

if any(idx<1)
    error('Elements in the second argument must be positive!')
end

i = zeros(n, 1);
i(idx) = 1;