function g = ahrens_dieter(a, b)

% Returns gamma variates, see Devroye (1986) page 410.
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
bb = (exp(1)+a)/exp(1);
cc = 1./a;
INDEX = 1:mm;
index = INDEX;
UW = NaN(nn,2);
V  = NaN(nn,1);
X  = NaN(nn,1);

while mm
    UW(index,:) = rand(mm,2);
    V(index) = UW(index,1).*bb(index);
    state1 = find(V(index)<=1);
    state2 = find(V(index)>1);
    ID = [];
    if ~isempty(state1)
        X(index(state1)) = V(index(state1)).^cc(index(state1));
        ID = INDEX(index(state1(UW(index(state1),2)>exp(-X(index(state1))))));
    end
    if ~isempty(state2)
        X(index(state2)) = -log(cc(index(state2)).*(bb(index(state2))-V(index(state2))));
        if isempty(ID)
            ID = INDEX(index(state2(UW(index(state2),2)>X(index(state2)).^(a(index(state2))-1))));
        else
            ID = [ID, INDEX(index(state2(UW(index(state2),2)>X(index(state2)).^(a(index(state2))-1))))];
        end
    end
    mm = length(ID);
    if mm
        index = ID;
    end
end

g = X.*b ;
