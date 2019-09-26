function rval = strjoin (cstr, delimiter)

% Adapted from Octave's implementation of strjoin
%
% Limitation: escaped characters (e.g. '\n') in delimiters will not be
% interpreted as the characters they represent.

% Copyright (C) 2013-2019 Ben Abbott
% Copyright (C) 2007 Muthiah Annamalai
% Copyright (C) 2019 Dynare Team
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

if nargin == 1
    delimiter = ' ';
elseif nargin < 1 || nargin > 2
    error('strjoin: must have either one or two arguments')
end
if ~(iscellstr(cstr) && (ischar(delimiter) || iscellstr(delimiter)))
    error('strjoin: first argument must be a cell array, second array either a char array or a cell array')
end

if numel(cstr) == 1
    rval = cstr{1};
    return;
end

if ischar(delimiter)
                        % There is no equivalent to do_string_escapes in MATLAB
                        %delimiter = do_string_escapes(delimiter);
    delimiter = {delimiter};
end

num = numel(cstr);
if numel(delimiter) == 1 && num > 1
    delimiter = repmat(delimiter, 1, num);
    delimiter(end) = {''};
elseif num > 0 && numel(delimiter) ~= num - 1
    error('strjoin: the number of delimiters does not match the number of strings');
else
    delimiter(end+1) = {''};
end

if num == 0
    rval = '';
else
    tmp = [cstr(:).'; delimiter(:).'];
    rval = [tmp{:}];
end

