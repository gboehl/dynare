function [fval, fjac] = broydentridiagonal(x)

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
    for k=1:n
        tmp = (3-2*x(k))*x(k);
        tmp1 = 0;
        if k~=1
            tmp1 = x(k-1);
        end
        tmp2 = 0;
        if k~=n
            tmp2 = x(k+1);
        end
        fval(k) = tmp-tmp1-2*tmp2+1;
    end
    if nargout>1
        fjac = zeros(n);
        for k=1:n
            fjac(k,k) = 3-4*x(k);
            if k~=1
                fjac(k,k-1) = -1;
            end
            if k~=n
                fjac(k,k+1) = -2;
            end
        end
    end
elseif isnumeric(x) && isvector(x) && all(isnan(x))
    if nargout==1
        fval = -ones(n, 1);
    else
        error('One output is required for initialization mode.')
    end
else
    error('Wrong input argument.')
end