function tf = contains(string, pattern, varargin)

% CONTAINS Returns 1 if the pattern is found in string, and 0 otherwise.
%
% INPUTS
% - string      [string, char, cell(str)]  String to be searhced.
% - pattern     [string, char, cell(str)]  The searched pattern.
%
% If 'IgnoreCase',IGNORE is provided, the function ignores case if IGNORE is true.
% The default value is false.
%
% OUTPUT
% - tf   [logical]
%
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
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

if mod(nargin,2) ~= 0
    error('contains: inputs should be a text, followed by the queried patterns, then Name-Value pair arguments.');
end

if ~((ischar(string) || isstring(string) || iscellstr(string)) && (ischar(pattern) || isstring(pattern) || iscellstr(pattern)))
    error('contains: first and second input arguments must be a string array, char array or cell array');
end

if ~isempty(varargin)
    if ~ischar(varargin{1})
        error('The #%d and #%d inputs must be a Name-Value pair.', nargin-1, nargin);
    end
    parameters = {'ignorecase', 'caseignore','insensitivecase', 'caseinsensitive','insensitive','ignore',...
                  'casesensitive', 'sensitive', 'sensitivecase'};
    if ~cellfun('isempty',(cellfun(@(s)strcmpi(varargin{1}, s), parameters, 'uni', 0)))
        case_ignore = varargin{end};
    else
        error('Unsupported parameter "%s".', varargin{1});
    end
else
    case_ignore = false;
end

string = cellstr(string);
pattern = cellstr(pattern);

if case_ignore
    string = lower(string);
    pattern = lower(pattern);
end

tf = false(size(string));
for ii = 1:numel(pattern)
    idx = regexp(string, pattern{ii});
    for jj = 1:numel(string)
        tf(jj) = tf(jj) || ~isempty(idx{jj});
    end
end
end