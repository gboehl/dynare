function  g = cheng(a, b)

% Returns gamma variates, see Devroye (1986) page 413.
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
bb = a-log(4);
cc = a+sqrt(2*a-1);
UV = NaN(nn,2);
Y  = NaN(nn,1);
X  = NaN(nn,1);
Z  = NaN(nn,1);
R  = NaN(nn,1);
index = 1:nn;
INDEX = index;

while mm
    UV(index,:) = rand(mm,2);
    Y(index) = a(index).*log(UV(index,2)./(1-UV(index,2)));
    X(index) = a(index).*exp(UV(index,2));
    Z(index) = UV(index,1).*UV(index,2).*UV(index,2);
    R(index) = bb(index) + cc(index).*Y(index)-X(index);
    jndex = index(R(index) >= 4.5*Z(index)-(1+log(4.5)));
    Jndex = setdiff(index, jndex);
    if ~isempty(Jndex)
        lndex = Jndex(R(Jndex) >= log(Z(Jndex)));
        Lndex = setdiff(Jndex, lndex);
    else
        Lndex = [];
    end
    index = INDEX(Lndex);
    mm = length(index);
end

g = X.*b;