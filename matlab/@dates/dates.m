function dd = dates(varargin) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function File} {@var{dd} =} dates (@var{a},@var{b},...)
%! @anchor{dates}
%! @sp 1
%! Constructor for the Dynare dates class (unordered sequence of dates).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item a
%! String, date.
%! @item b
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item dd
%! Dynare dates object.
%! @end table
%! @sp 1
%! @strong{Properties}
%! @sp 1
%! The constructor defines the following properties:
%! @sp 1
%! @table @ @var
%! @item ndate
%! Scalar integer, the number of dates.
%! @item freq
%! Scalar integer, the frequency of the time series. @var{freq} is equal to 1 if data are on a yearly basis or if
%! frequency is unspecified. @var{freq} is equal to 4 if data are on a quaterly basis. @var{freq} is equal to
%! 12 if data are on a monthly basis. @var{freq} is equal to 52 if data are on a weekly basis.
%! @item time
%! Array of integers (nobs*2). The first column defines the years associated to each date. The second column,
%! depending on the frequency, indicates the week, month or quarter numbers. For yearly data or unspecified frequency
%! the second column is filled by ones.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2011-2013 Dynare Team
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

if nargin>0 && ischar(varargin{1}) && isequal(varargin{1},'initialize')
    dd = struct('ndat', 0, 'freq', NaN(0), 'time', NaN(0,2));
    dd = class(dd,'dates');
    assignin('base','emptydatesobject',dd);
    return
end

dd = evalin('base','emptydatesobject');

switch nargin
  case 0
    % Returns an empty object
    return
  case 1
    if isa(varargin{1},'dates')
        % Returns a copy of the input argument
        dd = varargin{1};
    elseif isdate(varargin{1})
        date = string2date(varargin{1});
        dd.ndat = 1;
        dd.freq = date.freq;
        dd.time = date.time;
    elseif isfreq(varargin{1})
        % Instantiate an empty dates object (only set frequency)
        if ischar(varargin{1})
            dd.freq = string2freq(varargin{1});
        else
            dd.freq = varargin{1};
        end
    else
        error('dates:: Wrong calling sequence of the constructor!')
    end
  otherwise
    if isdate(varargin{1})
        dd.ndat = nargin;
        dd.time = NaN(dd.ndat,2);
        date = string2date(varargin{1});
        dd.freq = date.freq;
        dd.time(1,:) = date.time;
    elseif isdates(varargin{1})
        dd = horzcat(varargin{:});
        return
    elseif isfreq(varargin{1})
        S.type = '()';
        S.subs = varargin;
        dd = subsref(dd,S);
        return
    else
        error('dates::dates: Wrong calling sequence!')
    end
    for i=2:dd.ndat
        if isdate(varargin{i})
            date = string2date(varargin{i});
            if isequal(date.freq,dd.freq)
                dd.time(i,:) = date.time;
            else
                 error(sprintf('dates::dates: Check that all the inputs have the same frequency (see input number %i)!',i))
            end
        else
            error(sprintf('dates::dates: Input %i has to be a string date!',i))
        end
    end
end

%@test:1
%$ % Define some dates
%$ B1 = '1945Q3';
%$ B2 = '1950Q2';
%$ B3 = '1950q1';
%$ B4 = '1953Q4';
%$
%$ % Define expected results.
%$ e.time = [1945 3; 1950 2; 1950 1; 1953 4];
%$ e.freq = 4;
%$ e.ndat = 4;
%$
%$ % Call the tested routine.
%$ d = dates(B1,B2,B3,B4);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(d.time,e.time);
%$ t(2) = dyn_assert(d.freq,e.freq);
%$ t(3) = dyn_assert(d.ndat,e.ndat);
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define some dates
%$ B1 = '1945M3';
%$ B2 = '1950M2';
%$ B3 = '1950M10';
%$ B4 = '1953M12';
%$
%$ % Define expected results.
%$ e.time = [1945 3; 1950 2; 1950 10; 1953 12];
%$ e.freq = 12;
%$ e.ndat = 4;
%$
%$ % Call the tested routine.
%$ d = dates(B1,B2,B3,B4);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(d.time,e.time);
%$ t(2) = dyn_assert(d.freq,e.freq);
%$ t(3) = dyn_assert(d.ndat,e.ndat);
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define some dates
%$ B1 = '1945y';
%$ B2 = '1950Y';
%$ B3 = '1950a';
%$ B4 = '1953A';
%$
%$ % Define expected results.
%$ e.time = [1945 1; 1950 1; 1950 1; 1953 1];
%$ e.freq = 1;
%$ e.ndat = 4;
%$
%$ % Call the tested routine.
%$ d = dates(B1,B2,B3,B4);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(d.time,e.time);
%$ t(2) = dyn_assert(d.freq,e.freq);
%$ t(3) = dyn_assert(d.ndat,e.ndat);
%$ T = all(t);
%@eof:3

%@test:4
%$ % Define a dates object
%$ B = dates('1950Q1'):dates('1960Q3');
%$
%$
%$ % Call the tested routine.
%$ d = B(2);
%$ if isa(d,'dates')
%$     t(1) = 1;
%$ else
%$     t(1) = 0;
%$ end
%$
%$ if t(1)
%$     t(2) = dyn_assert(d.freq,B.freq);
%$     t(3) = dyn_assert(d.time,[1950 2]);
%$ end
%$ T = all(t);
%@eof:4

%@test:5
%$ % Define a dates object
%$ B = dates(4,1950,1):dates(4,1960,3);
%$
%$ % Call the tested routine.
%$ d = B(2);
%$ if isa(d,'dates')
%$     t(1) = 1;
%$ else
%$     t(1) = 0;
%$ end
%$
%$ if t(1)
%$     t(2) = dyn_assert(d.freq,B.freq);
%$     t(3) = dyn_assert(d.time,[1950 2]);
%$ end
%$ T = all(t);
%@eof:5

%@test:6
%$ % Define a dates object
%$ B = dates(4,[1950 1]):dates(4,[1960 3]);
%$
%$ % Call the tested routine.
%$ d = B(2);
%$ if isa(d,'dates')
%$     t(1) = 1;
%$ else
%$     t(1) = 0;
%$ end
%$
%$ if t(1)
%$     t(2) = dyn_assert(d.freq,B.freq);
%$     t(3) = dyn_assert(d.time,[1950 2]);
%$ end
%$ T = all(t);
%@eof:6