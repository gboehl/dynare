function CDF = chi2cdf(x, n)
% chi2cdf  CDF of the Chi2 distribution
%  CDF = chi2cdf(x, n) computes, for each element of X, the
%  CDF at X of the chi-square distribution with N degrees of freedom.
% Original file: statistics/distributions/chi2inv.m

% Copyright © 2013 Dynare Team
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

if (nargin ~= 2)
    error ('chi2cdf: you must give two arguments');
end

if (~isscalar (n))
    [retval, x, n] = common_size (x, n);
    if (retval > 0)
        error ('chi2cdf: x and n must be of common size or scalar');
    end
end

CDF = gamcdf (x, n / 2, 2);

end
