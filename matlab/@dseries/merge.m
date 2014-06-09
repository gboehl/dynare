function A = merge(B,C) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function File} {@var{A} =} merge (@var{B},@var{C})
%! @anchor{@dseries/plus}
%! @sp 1
%! Overloads the merge method for the Dynare time series class (@ref{dseries}).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item B
%! Dynare time series object instantiated by @ref{dseries}.
%! @item C
%! Dynare time series object instantiated by @ref{dseries}.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item A
%! Dynare time series object.
%! @end deftypefn
%@eod:

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

if ~isdseries(C)
    error('dseries::merge: Both inputs must be dseries objects!')
end

if ~isequal(frequency(B),frequency(C))
    error(['dseries::merge: Cannot merge ' inputname(1) ' and ' inputname(2) ' (frequencies are different)!'])
end

A = dseries();
[A.name, IBC, junk] = unique([B.name; C.name], 'last');
tex = [B.tex; C.tex];
A.tex = tex(IBC); 
A.vobs=length(IBC);

if B.nobs == 0
    A = C;
elseif C.nobs == 0
    A = B;
elseif firstdate(B) >= firstdate(C)
    diff = firstdate(B) - firstdate(C);
    A.nobs = max(B.nobs + diff, C.nobs);
    A.data = NaN(A.nobs, A.vobs);
    Z1 = [NaN(diff,B.vobs);B.data];
    if A.nobs > B.nobs + diff
        Z1 = [Z1; NaN(A.nobs-(B.nobs + diff),B.vobs)];
    end;
    Z2 = C.data;
    if A.nobs > C.nobs
        Z2 = [Z2; NaN(A.nobs - C.nobs,C.vobs)];
    end;
    Z = [Z1 Z2];
    A.data = Z(:,IBC);
    A_init = firstdate(C);
else
    diff = firstdate(C) - firstdate(B);
    A.nobs = max(C.nobs + diff, B.nobs);
    A.data = NaN(A.nobs, A.vobs);
    Z1 = [NaN(diff,C.vobs);C.data];
    if A.nobs > C.nobs + diff
        Z1 = [Z1; NaN(A.nobs-(C.nobs + diff),C.vobs)];
    end;
    Z2 = B.data;
    if A.nobs > B.nobs 
        Z2 = [Z2; NaN(A.nobs - B.nobs,B.vobs)];
    end;
    Z = [Z2 Z1];
    A.data = Z(:,IBC);
    A_init = B.init;
end

A.dates = A_init:A_init+(A.nobs-1);

%@test:1
%$ % Define a datasets.
%$ A = rand(10,2); B = randn(10,1);
%$
%$ % Define names
%$ A_name = {'A1';'A2'}; B_name = {'A1'};
%$
%$ t = zeros(4,1);
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dseries(A,[],A_name,[]);
%$    ts2 = dseries(B,[],B_name,[]);
%$    ts3 = merge(ts1,ts2);
%$    t(1) = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if length(t)>1
%$    t(2) = dyn_assert(ts3.vobs,2);
%$    t(3) = dyn_assert(ts3.nobs,10);
%$    t(4) = dyn_assert(ts3.data,[B, A(:,2)],1e-15);
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define a datasets.
%$ A = rand(10,2); B = randn(10,1);
%$
%$ % Define names
%$ A_name = {'A1';'A2'}; B_name = {'B1'};
%$
%$ t = zeros(4,1);
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dseries(A,[],A_name,[]);
%$    ts2 = dseries(B,[],B_name,[]);
%$    ts3 = merge(ts1,ts2);
%$    t(1) = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if length(t)>1
%$    t(2) = dyn_assert(ts3.vobs,3);
%$    t(3) = dyn_assert(ts3.nobs,10);
%$    t(4) = dyn_assert(ts3.data,[A, B],1e-15);
%$ end
%$ T = all(t);
%@eof:2