function [fval, fjac] = powell2(x)

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
    fval(1) = 10000*x(1)*x(2)-1;
    fval(2) = exp(-x(1))+exp(-x(2))-1.0001;
    if nargout>1
        fjac = zeros(2);
        fjac(1,1) = 10000*x(2);
        fjac(1,2) = 10000*x(1);
        fjac(2,1) = -exp(x(1));
        fjac(2,2) = -exp(x(2));
    end
elseif nargin
    if nargout==1
        fval = [0; 1];
    else
        error('One output is required for initialization mode.')
    end
else
    error('Wrong number of input arguments.')
end