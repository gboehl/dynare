function indices = kitagawa(weights, noise)

% Return indices for resampling.
%
% INPUTS
% - weights   [double]    n×1 vector of partcles' weights.
% - noise     [double]    scalar, uniform random deviates in [0,1]
%
% OUTPUTS
% - indices   [integer]   n×1 vector of indices in [1:n]

% Copyright © 2022-2023 Dynare Team
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

n= length(weights);

if nargin<2, noise = rand; end

indices = NaN(n, 1);

cweights = cumsum(weights);

wweights = (transpose(0:n-1)+noise)*(1.0/n);

j = 1;
for i=1:n
    while wweights(i)>cweights(j)
        j = j+1;
    end
    indices(i) = j;
end
