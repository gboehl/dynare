function  g = knuth(a, b)

% Returns gamma variates, see Bauwens, Lubrano & Richard (1999) page 316.
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
bb = sqrt(2*a-1);
dd = 1./(a-1);
Y = NaN(nn,1);
X = NaN(nn,1);
INDEX = 1:mm;
index = INDEX;

while mm
    Y(index) = tan(pi*rand(mm,1));
    X(index) = Y(index).*bb(index) + a(index) - 1 ;
    idy1 = find(X(index)>=0);
    idn1 = setdiff(index,index(idy1));
    if ~isempty(idy1)
        test = log(rand(length(idy1),1)) <= ...
               log(1+Y(index(idy1)).*Y(index(idy1))) + ...
               (a(index(idy1))-1).*log(X(index(idy1)).*dd(index(idy1))) - ...
               Y(index(idy1)).*bb(index(idy1)) ;
        idy2 = find(test);
        idn2 = setdiff(idy1, idy1(idy2));
    else
        idy2 = [];
        idn2 = [];
    end
    index = [ INDEX(idn1) , INDEX(index(idn2)) ] ;
    mm = length(index);
end

g = X.*b;