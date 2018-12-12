function g = weibull_rejection(a, b)

% Returns gamma variates, see Devroye (1986) page 415.
%
% INPUTS
% - a    [double]     n*1 vector, first hyperparameter.
% - b    [double]     n*1 vector, second hyperparameter.
%
% OUTPUTS
% - g    [double]     n*1 vector, gamma variates.

% Copyright (C) 2006-2018 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

nn = length(a);
mm = nn;
cc = 1./a ;
dd = a.^(a./(1-a)).*(1-a);
ZE = NaN(nn,2);
X  = NaN(nn,1);
INDEX = 1:mm;
index = INDEX;

while mm
    ZE(index,:) = exprnd(ones(mm,2));
    X(index) = ZE(index,1).^cc(index);
    id = find( (ZE(:,1)+ZE(:,2) > dd + X) );
    if isempty(id)
        mm = 0;
    else
        index = INDEX(id);
        mm = length(index);
    end
end

g = X.*b;