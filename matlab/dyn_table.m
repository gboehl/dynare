function dyn_table(title, preamble, afterward, rows, cols, indent, data)
%function dyn_table(title, rows, cols, indent, data)
% Print Table
%
% INPUTS
%   title      [char]
%   preamble   [cell string]
%   afterward  [cell string]
%   rows       [cell string]
%   cols       [cell string]
%   indent     [integer]
%   data       [matrix]
%
% OUTPUTS
%   None
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2017 Dynare Team
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

assert(ischar(title), 'title must be a char')
assert(iscellstr(preamble) && iscellstr(afterward) && iscellstr(rows) && iscellstr(cols), ...
    'preamble, afterward, rows, and cols must be cell arrays of strings')
assert(size(data, 1) == length(rows), 'must have the same number of rows')
assert(size(data, 2) == length(cols), 'must have the same number of columns')
assert(isint(indent), 'indent must be an integer')

skipline(3)

%% Print Output
rowstrlens = cellfun(@length, rows);
colstrlens = cellfun(@length, cols);
maxrowstrlen = max(rowstrlens);

colbegin = repmat(' ', 1, 2*indent + maxrowstrlen);
colrow = sprintf('%s', colbegin);
colrow2 = colrow;
format = ['    %-' num2str(maxrowstrlen) 's'];
for i = 1:length(cols)
    precision = 12;
    if colstrlens(i) < precision
        colrow = [colrow repmat(' ', 1, floor((precision-mod(colstrlens(i), precision))/2)) cols{i} repmat(' ', 1, ceil((precision-mod(colstrlens(i), precision))/2))];
        colrow2 = [colrow2 repmat('_', 1, precision)];
    else
        colrow = [colrow cols{i}];
        colrow2 = [colrow2 repmat('_', 1, colstrlens(i))];
        precision = colstrlens(i);
    end
    if i ~= length(cols)
        colrow = [colrow repmat(' ', 1, indent)];
        colrow2 = [colrow2 repmat(' ', 1, indent)];
    end
    format = [format repmat(' ', 1, indent) '%' sprintf('%d.5f', precision)];
end

% Center title
if length(title) >= length(colrow)
    fprintf('%s\n\n', title)
else
    fprintf('%s%s\n\n', repmat(' ', 1, floor((length(colrow)+indent-length(title))/2)), title);
end
spaces = repmat(' ', 1, indent);
for i = 1:length(preamble)
    fprintf('%s%s\n', spaces, preamble{i});
end

fprintf('%s\n', colrow);
fprintf('%s\n\n', colrow2);

format = [format '\n'];
for i = 1:length(rows)
    fprintf(format, rows{i}, data(i, :));
end

fprintf('\n');
for i = 1:length(afterward)
    fprintf('%s%s\n', spaces, afterward{i});
end

fprintf('%s\n\n', repmat('_', 1, length(colrow2)));
