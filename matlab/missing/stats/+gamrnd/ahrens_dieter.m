function g = ahrens_dieter(a, b)

% Returns gamma variates, see Devroye (1986) page 425.
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
bb = (exp(1)+a)/exp(1);
cc = 1./a;
index = 1:nn;
U = NaN(nn,1);
W = NaN(nn,1);
V = NaN(nn,1);
X = NaN(nn,1);

while mm
    U(index) = rand(mm,1);
    W(index) = rand(mm,1);
    V(index) = U(index).*bb(index);
    id1 = index(V(index)<=1);
    id2 = setdiff(index, id1);
    X(id1) = V(id1).^cc(id1);
    id3 = id1(W(id1)>exp(-X(id1)));
    X(id2) = -log(cc(id2).*(bb(id2)-V(id2)));
    id4 = id2(W(id2)>X(id2).^(a(id2)-1));
    index = [id3, id4];
    mm = length(index);
end

g = X.*b ;
