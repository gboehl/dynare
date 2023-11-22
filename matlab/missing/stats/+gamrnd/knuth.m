function  g = knuth(a, b)

% Returns gamma variates, see Knuth (1981) page 129.
%
% INPUTS
% - a    [double]     n*1 vector, first hyperparameter.
% - b    [double]     n*1 vector, second hyperparameter.
%
% OUTPUTS
% - g    [double]     n*1 vector, gamma variates.

% Copyright © 2006-2023 Dynare Team
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
bb = sqrt(2*a-1);
Y = NaN(nn,1);
X = NaN(nn,1);
index = 1:mm;

while mm
    Y(index) = tan(pi*rand(mm,1));
    X(index) = Y(index).*bb(index) + a(index) - 1;
    id1 = index(X(index)<=0); % Rejected draws.
    id2 = setdiff(index, id1);
    id3 = id2(rand(length(id2), 1)>(1+Y(id2).*Y(id2)).*exp((a(id2)-1).*(log(X(id2))-log(a(id2)-1))-bb(id2).*Y(id2))); % Rejected draws.
    index = [id1, id3];
    mm = length(index);
end

g = X.*b;
