function [l,c,d] = isint(a)
%  This function tests if the input argument is an integer.
%
%  INPUT
%    a    [double]   m*n matrix.
%
%  OUTPUT
%    b    [integer]  m*n matrix of 0 and 1. b(i,j)=1 if a(i,j) is an integer.
%    c    [integer]  p*1 vector of indices pointing to the integer elements of a.
%    d    [integer]  q*1 vector of indices pointing to the non integer elements of a.
%
%  SPECIAL REQUIREMENTS
%    None.
%
%  NOTES
%    p+q is equal to the product of m by n.

% Copyright © 2009-2017 Dynare Team
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

if ~isnumeric(a)
    l = false;
    if nargout>1
        c = [];
        d = [];
    end
    return
end

l = abs(fix(a)-a)<1e-15;

if nargout>1
    c = find(l==true);
    d = find(l==false);
end