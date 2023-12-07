function s=skewness(y)
% s=skewness(y)
% Compute normalized skewness of y
% Inputs:
%  - y  [double]  input vector
% Outputs:
%  - s  [double]  standardized skewness 

% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu

% Copyright © 2012 European Commission
% Copyright © 2012-2023 Dynare Team
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

m2=mean((y-mean(y)).^2);
m3=mean((y-mean(y)).^3);
s=m3/m2^1.5;