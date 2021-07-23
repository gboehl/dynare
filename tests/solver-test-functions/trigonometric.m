function [fval, fjac] = trigonometric(x)

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
        fval(j) = cos(x(j));
        sum = sum+fval(j);
    end
    for k=1:n
        fval(k) = (n+k)-sin(x(k))-sum-k*fval(k);
    end
    if nargout>1
        fjac = zeros(n);
        for j=1:n
            tmp = sin(x(j));
            for k=1:n
                fjac(k,j) = tmp;
            end
            fjac(j,j) = (j+1)*tmp-cos(x(j));
        end
    end
elseif isnumeric(x) && isvector(x) && all(isnan(x))
    if nargout==1
        fval = (1/n)*ones(n, 1);
    else
        error('One output is required for initialization mode.')
    end
else
    error('Wrong input argument.')
end