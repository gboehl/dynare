function [fval, fjac] = broydenbanded(x)

% Copyright © 2021 Dynare Team
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

n = length(x);

if isnumeric(x) && isvector(x) && all(~isnan(x))
    fval = zeros(n, 1);
    ml = 5;
    mu = 1;
    for k=1:n
        k1 = max(1, k-ml);
        k2 = min(k+mu, n);
        tmp = 0;
        for j=k1:k2
            if j~=k
                tmp = tmp + x(j)*(1+x(j));
            end
        end
        fval(k) = x(k)*(2+5*x(k)*x(k))+1-tmp;
    end
    if nargout>1
        fjac = zeros(n);
        ml = 5;
        mu = 1;
        for k=1:n
            k1 = max(1,k-ml);
            k2 = min(k+mu,n);
            for j=k1:k2
                if j~=k
                    fjac(k,j) = -(1+2*x(j));
                end
            end
            fjac(k,k) = 2+15*x(k)*x(k);
        end
    end
elseif isnumeric(x) && isvector(x) && all(isnan(x))
    if nargout==1
        fval = broydentridiagonal(nan(n,1));
    else
        error('One output is required for initialization mode.')
    end
else
    error('Wrong input argument.')
end