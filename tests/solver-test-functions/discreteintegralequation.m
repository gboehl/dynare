function [fval, fjac] = discreteintegralequation(x)

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
    h = 1.0/(n+1);
    for k=1:n
        tk = k*h;
        sum1 = 0;
        for j=1:k
            tj = j*h;
            tmp = (x(j)+tj+1)^3;
            sum1 = sum1+tj*tmp;
        end
        sum2 = 0;
        kp1 = k+1;
        if n>=kp1
            for j=kp1:n
                tj = j*h;
                tmp = (x(j)+tj+1)^3;
                sum2 = sum2+(1-tj)*tmp;
            end
        end
        fval(k) = x(k) + h*((1-tk)*sum1+tk*sum2)/2;
    end
    if nargout>1
        fjac = zeros(n);
        h = 1.0/(n+1);
        for k=1:n
            tk = k*h;
            for j=1:n
                tj = j*h;
                tmp = 3*(x(j)+tj+1)^2;
                fjac(k,j) = h*min(tj*(1-tk), tk*(1-tj))*tmp/2;
            end
            fjac(k,k) = fjac(k,k)+1;
        end
    end
elseif isnumeric(x) && isvector(x) && all(isnan(x))
    if nargout==1
        fval = discreteboundaryvalue(nan(n,1));
    else
        error('One output is required for initialization mode.')
    end
else
    error('Wrong input argument.')
end