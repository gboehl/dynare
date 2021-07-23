function [fval, fjac] = variablydimensioned(x)

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
    sum = 0;
    for j=1:n
        sum = sum+j*(x(j)-1);
    end
    tmp = sum*(1+2*sum*sum);
    for k=1:n
        fval(k) = x(k)-1+k*tmp;
    end
    if nargout>1
        fjac = zeros(n);
        sum = 0;
        for j=1:n
            sum = sum + j*(x(j)-1);
        end
        tmp = 1+6*sum*sum;
        for k=1:n
            for j=1:n
                fjac(k,j) = k*j*tmp;
                fjac(j,k) = fjac(k,j);
            end
            fjac(k,k) = fjac(k,k)+1;
        end
    end
elseif isnumeric(x) && isvector(x) && all(isnan(x))
    if nargout==1
        h = 1.0/n;
        fval = zeros(n, 1);
        for j=1:n
            fval(j) = 1-j*h;
        end
    else
        error('One output is required for initialization mode.')
    end
else
    error('Wrong input argument.')
end