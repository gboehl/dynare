function d=hess_element(func,arg1,arg2,elem1,elem2,args)
% function d=hess_element(func,arg1,arg2,elem1,elem2,args)
% returns an entry of the finite differences approximation to the hessian of func
%
% INPUTS
%    func       [function name]    string with name of the function
%    arg1       [int]              the indices showing the element within the hessian that should be returned
%    arg2       [int]
%    elem1      [int]              vector index 1
%    elem2      [int]              vector index 2
%    args       [cell array]       arguments provided to func
%
% OUTPUTS
%    d          [double]           the (arg1,arg2) entry of the hessian
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2010-2020 Dynare Team
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

assert(arg1 <= length(args) && arg2 <= length(args));

func = str2func(func);

h=1e-6;
p10 = args;
p01 = args;
m10 = args;
m01 = args;
p11 = args;
m11 = args;

p10{arg1}(elem1) = p10{arg1}(elem1) + h;
m10{arg1}(elem1) = m10{arg1}(elem1) - h;

p11{arg1}(elem1) = p11{arg1}(elem1) + h;
m11{arg1}(elem1) = m11{arg1}(elem1) - h;

p01{arg2}(elem2) = p01{arg2}(elem2) + h;
m01{arg2}(elem2) = m01{arg2}(elem2) - h;

p11{arg2}(elem2) = p11{arg2}(elem2) + h;
m11{arg2}(elem2) = m11{arg2}(elem2) - h;

% From Abramowitz and Stegun. Handbook of Mathematical Functions (1965)
% formulas 25.3.24 and 25.3.27 p. 884
if arg1==arg2
    d = (16*func(p10{:})...
         +16*func(m10{:})...
         -30*func(args{:})...
         -func(p11{:})...
         -func(m11{:}))/(12*h^2);
else
    d = (func(p10{:})...
         +func(m10{:})...
         +func(p01{:})...
         +func(m01{:})...
         -2*func(args{:})...
         -func(p11{:})...
         -func(m11{:}))/(-2*h^2);
end
end
