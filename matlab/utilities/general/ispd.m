function [test, penalty] = ispd(A)

%@info:
%! @deftypefn {Function File} {[@var{test}, @var{penalty}  =} ispd (@var{A})
%! @anchor{ispd}
%! @sp 1
%! Tests if the square matrix @var{A} is positive definite.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item A
%! A square matrix.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item test
%! Integer scalar equal to 1 if @var{A} is a positive definite sqquare matrix, 0 otherwise.
%! @item penalty
%! Absolute value of the uum of the negative eigenvalues of A. This output argument is optional.
%! @end table
%! @end deftypefn
%@eod:

% Copyright © 2007-2022 Dynare Team
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

if ~isquare(A)
    error(['ispd:: Input argument ' inputname(1) ' has to be a square matrix!'])
end

[~, info] = chol(A);
test = ~info;

if nargout>1
    penalty = 0;
    if info
        if isoctave && any(any(~isfinite(A))) % workaround for https://savannah.gnu.org/bugs/index.php?63082
            penalty = 1;
        else
            a = diag(eig(A));
            k = find(a<0);
            if k > 0
                penalty = sum(-a(k));
            end
        end
    end
end