function [fval, fjac] = helicalvalley(x)

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

C7 = 0.25;
C8 = 0.5;
TPI = 8*atan(1);

if nargin==1
    fval = zeros(3, 1);
    tmp1 = C7*sign(x(2));
    if x(1)>0
        tmp1 = atan(x(2)/x(1))/TPI;
    end
    if x(1)<0
        tmp1 = atan(x(2)/x(1))/TPI+C8;
    end
    tmp2 = sqrt(x(1)*x(1)+x(2)*x(2));
    fval(1) = 10*(x(3)-10*tmp1);
    fval(2) = 10*(tmp2-1);
    fval(3) = x(3);
    if nargout>1
        fjac = zeros(3);
        tmp = x(1)*x(1)+x(2)*x(2);
        tmp1 = TPI*tmp;
        tmp2 = sqrt(tmp);
        fjac(1,1) = 100*x(2)/tmp1;
        fjac(1,2) = -100*x(1)/tmp1;
        fjac(1,3) = 10;
        fjac(2,1) = 10*x(1)/tmp2;
        fjac(2,2) = 10*x(2)/tmp2;
        fjac(2,3) = 0;
        fjac(3,1) = 0;
        fjac(3,2) = 0;
        fjac(3,3) = 1;
    end
elseif ~nargin
    if nargout==1
        fval = [-1; 0; 0];
    else
        error('One output is required for initialization mode.')
    end
else
    error('Wrong number of input arguments.')
end