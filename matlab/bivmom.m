function [y,dy] = bivmom(p,rho)
% Computes the product moment (and its derivative with respect to standard
% errors and correlation parameters) of X_1^{p_1}X_2^{p_2}, where X_1 and X_2
% are standard bivariate normally distributed.
% n : dimension of X
% rho: correlation coefficient between X_1 and X_2
% =========================================================================
% INPUTS
%   p   [2 by 1]    powers of X_{1} and X_{2}
%   rho [1 by 1]    correlation coefficient between X_1 and X_2
% -------------------------------------------------------------------------
% OUTPUTS
%   y   [1 by 1]    product moment E[X_1^{p_1}X_2^{p_2}]
%   dy  [1 by 1]    derivative of y wrt to rho
% -------------------------------------------------------------------------
% This function is based upon bivmom.m which is part of replication codes
% of the following paper:
% Kan, R.: "From moments of sum to moments of product." Journal of 
% Multivariate Analysis, 2008, vol. 99, issue 3, pages 542-554.
% bivmom.m can be retrieved from http://www-2.rotman.utoronto.ca/~kan/papers/prodmom.zip
% Further references:
% Kotz, Balakrishnan, and Johnson (2000), Continuous Multivariate Distributions, Vol. 1, p.261
% Note that there is a typo in Eq.(46.25), there should be an extra rho in front 
% of the equation.
% =========================================================================
% Copyright (C) 2008-2015 Raymond Kan <kan@chass.utoronto.ca>
% Copyright (C) 2019-2020 Dynare Team
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
% =========================================================================
s1 = p(1);
s2 = p(2);
rho2 = rho^2;
if nargout > 1
    drho2 = 2*rho;
end
if rem(s1+s2,2)==1
   y = 0;
   return
end
r = fix(s1/2);
s = fix(s2/2);
y = 1;
c = 1;
if nargout > 1
    dy = 0;
    dc = 0;
end
odd = 2*rem(s1,2);
for j=1:min(r,s)
    if nargout > 1
        dc = 2*dc*(r+1-j)*(s+1-j)*rho2/(j*(2*j-1+odd)) + 2*c*(r+1-j)*(s+1-j)*drho2/(j*(2*j-1+odd));    
    end
    c = 2*c*(r+1-j)*(s+1-j)*rho2/(j*(2*j-1+odd));
    y = y+c;
    if nargout > 1
        dy = dy + dc;
    end
end
if odd
   if nargout > 1
    dy = y + dy*rho;
   end
   y = y*rho;
end
y = prod([1:2:s1])*prod([1:2:s2])*y;
if nargout > 1
    dy = prod([1:2:s1])*prod([1:2:s2])*dy;
end

