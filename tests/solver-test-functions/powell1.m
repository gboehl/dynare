function [fval, fjac] = powell1(x)

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

if nargin==1
    fval = zeros(4, 1);
    fval(1) = x(1)+10*x(2);
    fval(2) = sqrt(5)*(x(3)-x(4));
    tmp = x(2)-2*x(3);
    fval(3) = tmp*tmp;
    tmp = x(1)-x(4);
    fval(4) = sqrt(10)*tmp*tmp;
    if nargout>1
        fjac = zeros(4);
        fjac(1,1) = 1;
        fjac(1,2) = 10;
        fjac(2,3) = sqrt(5);
        fjac(2,4) = -fjac(2,3);
        fjac(3,2) = 2*(x(2)-2*x(3));
        fjac(3,3) = -2*fjac(3,2);
        fjac(4,1) = 2*sqrt(10)*(x(1)-x(4));
        fjac(4,4) = -fjac(4,1);
    end
elseif ~nargin
    if nargout==1
        fval = [3; -1; 0; 1];
    else
        error('One output is required for initialization mode.')
    end
else
    error('Wrong number of input arguments.')
end