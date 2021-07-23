function [fval, fjac] = watson(x)

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

C9 = 29.0;

if isnumeric(x) && isvector(x) && all(~isnan(x))
    fval = zeros(n, 1);
    for i=1:29
        ti = i/C9;
        sum1 = 0;
        tmp = 1;
        for j=2:n
            sum1 = sum1+(j-1)*tmp*x(j);
            tmp = tmp*ti;
        end
        sum2 = 0;
        tmp = 1;
        for j=1:n
            sum2 = sum2+tmp*x(j);
            tmp = tmp*ti;
        end
        tmp1 = sum1-sum2*sum2-1;
        tmp2 = 2*ti*sum2;
        tmp = 1/ti;
        for k=1:n
            fval(k) = fval(k) + tmp*((k-1)-tmp2)*tmp1;
            tmp = tmp*ti;
        end
    end
    tmp = x(2)-x(1)*x(1)-1;
    fval(1) = fval(1)+ x(1)*(1-2*tmp);
    fval(2) = fval(2)+tmp;
    if nargout>1
        fjac = zeros(n);
        for i=1:29
            ti = i/C9;
            sum1 = 0;
            tmp = 1;
            for j=2:n
                sum1 = sum1 + (j-1)*tmp*x(j);
                tmp = tmp*ti;
            end
            sum2 = 0;
            tmp = 1;
            for j=1:n
                sum2 = sum2+tmp*x(j);
                tmp = tmp*ti;
            end
            tmp1 = 2*(sum1-sum2*sum2-1);
            tmp2 = 2*sum2;
            tmp = ti*ti;
            tk = 1;
            for k=1:n
                tj = tk;
                for j=k:n
                    fjac(k,j) = fjac(k,j)+tj*(((k-1)/ti-tmp2)*((j-1)/ti-tmp2)-tmp1);
                    tj = tj*ti;
                end
                tk = tk*tmp;
            end
        end
        fjac(1,1) = fjac(1,1)+6*x(1)*x(1)-2*x(2)+3;
        fjac(1,2) = fjac(1,2)-2*x(1);
        fjac(2,2) = fjac(2,2)+1;
        for k=1:n
            for j=k:n
                fjac(j,k) = fjac(k,j);
            end
        end
    end
elseif isnumeric(x) && isvector(x) && all(isnan(x))
    if nargout==1
        fval = zeros(length(x), 1);
    else
        error('One output is required for initialization mode.')
    end
else
    error('Wrong input argument.')
end