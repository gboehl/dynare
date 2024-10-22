function  [ldens,Dldens,D2ldens] = lpdfnorm(x,a,b)
% Evaluates the logged UNIVARIATE GAUSSIAN PDF at x.
%
% INPUTS
%    x     [double]  m*n matrix of locations,
%    a     [double]  m*n matrix or scalar, First GAUSSIAN distribution parameters (expectation)
%    b     [double]  m*n matrix or scalar, Second GAUSSIAN distribution parameters (standard deviation).
%
% OUTPUTS
%    ldens   [double]  m*n matrix of logged GAUSSIAN densities evaluated at x.
%    Dldens  [double]  m*n matrix of first derivatives of logged GAUSSIAN densities.
%    D2ldens [double]  m*n matrix of second derivatives of logged  matrix of GAUSSIAN densities.
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2003-2021 Dynare Team
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

if nargin<3, b=1; end
if nargin<2, a=0; end
ldens = -log(b) -.5*log(2*pi) - .5*((x-a)./b).*((x-a)./b) ;

if nargout >1
    Dldens =  - (1./b).*((x-a)./b) ;
end

if nargout == 3
    D2ldens =  - (1./b).^2 ;
end