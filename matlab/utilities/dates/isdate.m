function b = isdate(str)  % --*-- Unitary tests --*--

% Tests if the input string can be interpreted as a date.
%
% INPUTS 
%  o str     string.
%
% OUTPUTS 
%  o b       integer scalar, equal to 1 if str can be interpreted as a date or 0 otherwise.

% Copyright (C) 2013 Dynare Team
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

%if length(str)>1
    b = isquaterly(str) || isyearly(str) || ismonthly(str) || isweekly(str);
    %else
    %b = 0;
    %end

%@test:1
%$
%$ date_1 = '1950M2';
%$ date_2 = '1950m2';
%$ date_3 = '-1950m2';
%$ date_4 = '1950m52';
%$ date_5 = ' 1950';
%$ date_6 = '1950Y';
%$ date_7 = '-1950a';
%$ date_8 = '1950m ';
%$
%$ t(1) = dyn_assert(isdate(date_1),1);
%$ t(2) = dyn_assert(isdate(date_2),1);
%$ t(3) = dyn_assert(isdate(date_3),1);
%$ t(4) = dyn_assert(isdate(date_4),0);
%$ t(5) = dyn_assert(isdate(date_5),0);
%$ t(6) = dyn_assert(isdate(date_6),1);
%$ t(7) = dyn_assert(isdate(date_7),1);
%$ t(8) = dyn_assert(isdate(date_8),0);
%$ T = all(t);
%@eof:1