function [fval, fjac] = wood(x)

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

C3 = 200;
C4 = 20.2;
C5 = 19.8;
C6 = 180.0;

if nargin==1
    fval = zeros(4, 1);
    tmp1 = x(2)-x(1)*x(1);
    tmp2 = x(4)-x(3)*x(3);
    fval(1) = -C3*x(1)*tmp1-(1-x(1));
    fval(2) = C3*tmp1+C4*(x(2)-1.0)+C5*(x(4)-1);
    fval(3) = -C6*x(3)*tmp2-(1-x(3));
    fval(4) = C6*tmp2+C4*(x(4)-1)+C5*(x(2)-1);
    if nargout>1
        fjac = zeros(4);
        tmp1 = x(2)-3*x(1)*x(1);
        tmp2 = x(4)-3*x(3)*x(3);
        fjac(1,1) = -C3*tmp1+1;
        fjac(1,2) = -C3*x(1);
        fjac(2,1) = -2*C3*x(1);
        fjac(2,2) = C3+C4;
        fjac(2,4) = C5;
        fjac(3,3) = -C6*tmp2+1;
        fjac(3,4) = -C6*x(3);
        fjac(4,2) = C5;
        fjac(4,3) = -2*C6*x(3);
        fjac(4,4) = C6+C4;
    end
elseif ~nargin
    if nargout==1
        fval = [-3; -1; -3; -1];
    else
        error('One output is required for initialization mode.')
    end
else
    error('Wrong number of input arguments.')
end