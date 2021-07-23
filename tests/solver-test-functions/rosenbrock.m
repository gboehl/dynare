function [fval, fjac] = rosenbrock(x)

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
    fval = zeros(2, 1);
    fval(1) = 1.0-x(1);
    fval(2) = 10.0*(x(2)-fval(1)*fval(1));
    if nargout>1
        fjac = zeros(2);
        fjac(1,1) = -1;
        fjac(1,2) = 0;
        fjac(2,1) = -20*x(1);
        fjac(2,2) = 10;
    end
elseif nargin
    if nargout==1
        fval = [10000; 1];
    else
        error('One output is required for initialization mode.')
    end
else
    error('Wrong number of input arguments.')
end