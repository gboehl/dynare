function g = best_1983(a, b)

% Returns gamma variates, see Devroye (1986) page 426 and Best (1983) page 187.
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
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

nn = length(a);
mm = nn;
index = 1:nn;
U = NaN(nn,1);
Ustar = NaN(nn, 1);
P = NaN(nn,1);
X = NaN(nn,1);
Y = NaN(nn,1);
zz = .07 + .75*sqrt(1-a);
bb = 1 + exp(-zz).*a./zz;
cc = 1./a;

while mm
    U(index) = rand(mm,1);
    Ustar(index) = rand(mm, 1);
    P(index) = U(index).*bb(index);
    id1 = index(P(index)<=1);
    id2 = setdiff(index, id1); % Goto 4.
    X(id1) = zz(id1).*(P(id1).^cc(id1));
    id3 = id1(Ustar(id1)>((2-X(id1))./(2+X(id1))));
    id5 = id3(Ustar(id3)>exp(-X(id3)));
    X(id2) = -log(cc(id2).*zz(id2).*(bb(id2)-P(id2))); % This is 4.
    Y(id2) = X(id2)./zz(id2);
    id4 = id2(Ustar(id2).*(a(id2)+(1-a(id2)).*Y(id2))>1);
    id6 = id4(Ustar(id4)>Y(id4).^(a(id4)-1));
    index = [id5, id6];
    mm = length(index);
end

g = X.*b;
