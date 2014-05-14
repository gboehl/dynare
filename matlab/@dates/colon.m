function C = colon(varargin) % --*-- Unitary tests --*--

% Overloads the colon operator (:). This method can be used to create ranges of dates.
%
% INPUTS 
%  o A    dates object with one element.
%  o d    integer scalar, number of periods between each date (default value, if nargin==2, is one)
%  o B    dates object with one element.
%
% OUTPUTS 
%  o C    dates object with length(B-A) elements (if d==1).
%
% REMARKS 
%  B must be greater than A if d>0.

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

if isequal(nargin,2)    
    A = varargin{1};
    B = varargin{2};
    d = 1;
    if ~(isa(A,'dates') && isa(B,'dates') && isequal(length(A),1) && isequal(length(B),1))
        error('dates::colon: In an expression like A:B, A and B must be dates objects!')
    end
elseif isequal(nargin,3)
    A = varargin{1};
    B = varargin{3};
    d = varargin{2};
    if ~(isa(A,'dates') && isa(B,'dates') && isequal(length(A),1) && isequal(length(B),1))
        error('dates::colon: In an expression like A:d:B, A and B must be dates objects and d a scalar integer (number of periods)!')
    end
    if ~(isscalar(d) && isint(d))
        error('dates::colon: In an expression like A:d:B, A and B must be dates objects and d a scalar integer (number of periods)!')
    end
    if isequal(d,0)
        error('dates::colon: In an expression like A:d:B, d (the incremental number of periods) must nonzero!')
    end
else
    error('dates::colon: Wrong calling sequence! See the manual for the colon (:) operator and dates objects.')
end

if ~isequal(A.freq,B.freq)
    error(['dates::colon: Input arguments ' inputname(1) 'and ' inputname(2) ' must have common frequency!'])
end

if A>B && d>0
    error(['dates::colon: ' inputname(1) ' must precede ' inputname(2) '!' ])
end

if B>A && d<0
    error(['dates::colon: ' inputname(2) ' must precede ' inputname(1) '!' ])
end

C = dates();
n = (B-A)+1;
if d>1
    n = length(1:d:n);
end
C.freq = A.freq;
C.ndat = n;
C.time = NaN(n,2);
C.time(1,:) = A.time;

for linee=2:n
    C.time(linee,:) = add_periods_to_date(C.time(linee-1,:), C.freq, d) ;%add_periods_to_array_of_dates(C.time(linee-1,:), C.freq, d);
end

%@test:1
%$ % Define two dates
%$ date_1 = '1950Q2';
%$ date_2 = '1951Q4';
%$
%$ % Define expected results.
%$ e.freq = 4;
%$ e.time = [1950 2; 1950 3; 1950 4; 1951 1; 1951 2; 1951 3; 1951 4];
%$
%$ % Call the tested routine.
%$ d1 = dates(date_1);
%$ d2 = dates(date_2);
%$ d3 = d1:d2;
%$
%$ % Check the results.
%$ t(1) = dyn_assert(d3.time,e.time);
%$ t(2) = dyn_assert(d3.freq,e.freq);
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define expected results.
%$ e.freq = 4;
%$ e.time = [1950 2; 1950 3; 1950 4; 1951 1; 1951 2; 1951 3; 1951 4];
%$
%$ % Call the tested routine.
%$ d = dates('1950Q2'):dates('1951Q4');
%$
%$ % Check the results.
%$ t(1) = dyn_assert(d.time,e.time);
%$ t(2) = dyn_assert(d.freq,e.freq);
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define expected results.
%$ e.freq = 4;
%$ e.time = [1950 2; 1950 4; 1951 2; 1951 4];
%$
%$ % Call the tested routine.
%$ d = dates('1950Q2'):2:dates('1951Q4');
%$
%$ % Check the results.
%$ t(1) = dyn_assert(d.time,e.time);
%$ t(2) = dyn_assert(d.freq,e.freq);
%$ T = all(t);
%@eof:3


%$ @test:3
%$ % Create an empty dates object for quaterly data
%$ qq = dates('Q');
%$
%$ % Define expected results.
%$ e.freq = 4;
%$ e.time = [1950 2; 1950 3; 1950 4; 1951 1; 1951 2; 1951 3; 1951 4];
%$
%$ % Call the tested routine.
%$ d = qq(1950,2):qq(1951,4);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(d.time,e.time);
%$ t(2) = dyn_assert(d.freq,e.freq);
%$ T = all(t);
%$ @eof:3