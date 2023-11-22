function  g = best_1978(a ,b)

% Returns gamma variates, see Devroye (1986) page 410 and Best (1978).
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
bb = a-1;
cc = 3*a-.75;
U = NaN(nn,1);
Y = NaN(nn,1);
X = NaN(nn,1);
Z = NaN(nn,1);
W = NaN(nn,1);
index = 1:nn;

while mm
    U(index) = rand(mm,1);
    W(index) = U(index).*(1-U(index)); % e
    Y(index) = sqrt(cc(index)./W(index)).*(U(index)-.5); % f
    X(index) = bb(index)+Y(index); % x
    id1 = index(X(index)<0); % Reject.
    id2 = setdiff(index, id1);
    Z(id2) = 64.0*(W(id2).^3).*(rand(length(id2),1).^2); % d
    id3 = id2(Z(id2)>1.0-2.0*Y(id2).*Y(id2)./X(id2)); % Reject.
    id4 = id3(log(Z(id3))>2.0*(bb(id3).*log(X(id3)./bb(id3))-Y(id3))); % Reject.
    index = [id1, id4];
    mm = length(index);
end

g = X.*b;
