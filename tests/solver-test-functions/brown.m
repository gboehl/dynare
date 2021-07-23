function [fval, fjac] = brown(x)

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
    sum = -(n+1);
    prod = 1;
    for j=1:n
        sum = sum+x(j);
        prod = prod*x(j);
    end
    for k=1:n-1
        fval(k) = x(k)+sum;
    end
    fval(n) = prod-1;
    if nargout>1
        fjac = zeros(n);
        prod = 1;
        for j=1:n
            prod = prod*x(j);
            for k=1:n
                fjac(k,j) = 1;
            end
            fjac(j,j) = 2;
        end
        for j=1:n
            tmp = x(j);
            if abs(tmp)<eps()
                tmp = 1;
                prod = 1;
                for k=1:n
                    if k~=j
                        prod = prod*x(k);
                    end
                end
            end
            fjac(n,j) = prod/tmp;
        end
    end
elseif isnumeric(x) && isvector(x) && all(isnan(x))
    if nargout==1
        fval = 0.5*ones(length(x), 1); 
    else
        error('One output is required for initialization mode.')
    end
else
    error('Wrong input argument.')
end