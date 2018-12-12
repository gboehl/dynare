function  g = best_1978(a ,b)

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
bb = a-1;
cc = 3*a-.75;
UV = NaN(nn,2);
Y  = NaN(nn,1);
X  = NaN(nn,1);
Z  = NaN(nn,1);
W  = NaN(nn,1);
index = 1:nn;
INDEX = index;

while mm
    UV(index,:) = rand(mm,2);
    W(index) = UV(index,1).*(1-UV(index,1));
    Y(index) = sqrt(cc(index)./W(index)).*(UV(index,1)-.5);
    X(index) = bb(index)+Y(index);
    jndex = index(X(index)>=0);
    Jndex = setdiff(index,jndex);
    if ~isempty(jndex)
        Z(jndex) = 64*W(jndex).*W(jndex).*W(jndex).*UV(jndex,2).*UV(jndex,2);
        kndex = jndex(Z(jndex)<=1-2*Y(jndex).*Y(jndex)./X(jndex));
        Kndex = setdiff(jndex, kndex);
        if ~isempty(Kndex)
            lndex = Kndex(log(Z(Kndex))<=2*(bb(Kndex).*log(X(Kndex)./bb(Kndex))-Y(Kndex)));
            Lndex = setdiff(Kndex, lndex);
        else
            Lndex = [];
        end
        new_index = INDEX(Lndex);
    end
    index = union(new_index, INDEX(Jndex));
    mm = length(index);
end

g = X.*b;
