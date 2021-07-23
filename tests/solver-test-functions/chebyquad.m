function [fval, fjac] = chebyquad(x)

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
    for j=1:n
        tmp1 = 1;
        tmp2 = 2*x(j)-1;
        tmp = 2*tmp2;
        for i=1:n
            fval(i) = fval(i)+tmp2;
            ti = tmp*tmp2-tmp1;
            tmp1 = tmp2;
            tmp2 = ti;
        end
    end
    tk = 1.0/n;
    iev = -1;
    for k=1:n
        fval(k) = fval(k)*tk;
        if iev>0
            fval(k) = fval(k)+1/(k^2-1);
        end
        iev = -iev;
    end
    if nargout>1
        fjac = zeros(n);
        tk = 1.0/n;
        for j=1:n
            tmp1 = 1;
            tmp2 = 2*x(j)-1;
            tmp = 2*tmp2;
            tmp3 = 0;
            tmp4 = 2;
            for k=1:n
                fjac(k,j) = tk*tmp4;
                ti = 4*tmp2+tmp*tmp4-tmp3;
                tmp3 = tmp4;
                tmp4 = ti;
                ti = tmp*tmp2-tmp1;
                tmp1 = tmp2;
                tmp2 = ti;
            end
        end
    end
elseif isnumeric(x) && isvector(x) && all(isnan(x))
    if nargout==1
        h = 1.0/(n+1);
        fval = zeros(length(x), 1); 
        for j=1:n
            fval(j) = j*h;
        end
    else
        error('One output is required for initialization mode.')
    end
else
    error('Wrong input argument.')
end