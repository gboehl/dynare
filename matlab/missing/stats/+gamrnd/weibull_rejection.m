function g = weibull_rejection(a, b)

% Returns gamma variates, using Weibull rejection algorithm, see Väduva (1977) page 556.
%
% INPUTS
% - a    [double]     n*1 vector, first hyperparameter.
% - b    [double]     n*1 vector, second hyperparameter.
%
% OUTPUTS
% - g    [double]     n*1 vector, gamma variates.

% Copyright © 2006-2018 Dynare Team
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

nn = length(a);
mm = nn;
cc = 1./a ;
zz = a.^(a./(1-a));
aa = exp(zz.*(a-1));
Z = NaN(nn, 1);
Y = NaN(nn, 1);
X = NaN(nn, 1);
index = 1:nn;

while mm
    Z(index) = -log(rand(mm, 1));
    Y(index) = Z(index).^cc(index);
    INDEX = index(rand(mm,1)>aa(index).*exp(Z(index)-Y(index)));
    id = setdiff(index, INDEX);
    X(id) = Y(id);
    index = INDEX;
    mm = length(index);
end

g = X.*b;