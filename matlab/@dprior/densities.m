function lpd = densities(o, X)

% Evaluate the logged prior densities at X (for each column).
%
% INPUTS
% - o       [dprior]
% - X       [double]   m×n matrix, n points where the prior density is evaluated.
%
% OUTPUTS
% - lpd     [double]   1×n, values of the logged prior density at X.

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

n = columns(X);

lpd = NaN(1, n);

parfor i=1:n
    lpd(i) = density(o, X(:,i));
end
