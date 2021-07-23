function [fval, fjac] = discreteboundaryvalue(x)

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
    hh = h*h;
    for k=1:n
        tmp = (x(k)+k*h+1)^3;
        tmp1 = 0;
        if k~=1
            tmp1 = x(k-1);
        end
        tmp2 = 0;
        if k~=n
            tmp2 = x(k+1);
        end
        fval(k) = 2*x(k)-tmp1-tmp2+tmp*hh/2;
    end
    if nargout>1
        fjac = zeros(n);
        h = 1/(n+1);
        hh = h*h;
        for k=1:n
            tmp = 3*(x(k)+k*h+1)^2;
            fjac(k,k) = 2+tmp*hh/2;
            if k~=1
                fjac(k,k-1) = -1;
            end
            if k~=n
                fjac(k,k+1) = -1;
            end
        end
    end
elseif isnumeric(x) && isvector(x) && all(isnan(x))
    if nargout==1
        h = 1.0/(n+1);
        fval = zeros(n, 1);
        for j=1:n
            tj = j*h;
            fval(j) = tj*(tj-1);
        end
    else
        error('One output is required for initialization mode.')
    end
else
    error('Wrong input argument.')
end