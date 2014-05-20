function o = writeSeriesForTable(o, fid, dates, precision)
%function o = writeSeriesForTable(o, fid, dates, precision)
% Write Table Row
%
% INPUTS
%   o            [report_series]    report_series object
%   fid          [int]              file id
%   dates        [dates]            dates for report_series slice
%   precision    [float]            precision with which to print the data
%
%
% OUTPUTS
%   o            [report_series]    report_series object
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2013-2014 Dynare Team
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

%% Validate options passed to function
assert(fid ~= -1);
for i=1:length(dates)
    assert(isdates(dates{i}));
end
assert(isint(precision));

%% Validate options provided by user
assert(ischar(o.tableSubSectionHeader), '@report_series.writeSeriesForTable: tableSubSectionHeader must be a string');
if isempty(o.tableSubSectionHeader)
    assert(~isempty(o.data) && isdseries(o.data), ...
           '@report_series.writeSeriesForTable: must provide data as a dseries');

    if ~isempty(o.tableDataRhs)
        assert(~isempty(o.tableDataRhs) && isdseries(o.tableDataRhs), ...
               '@report_series.writeSeriesForTable: must provide tableDataRhs as a dseries');
        assert(iscell(dates) && length(dates) == 2, ...
               '@report_series.writeSeriesForTable: must provide second range with tableDataRhs');
    end
end

assert(ischar(o.tableNegColor), '@report_series.writeSeriesForTable: tableNegColor must be a string');
assert(ischar(o.tablePosColor), '@report_series.writeSeriesForTable: tablePosColor must be a string');
assert(ischar(o.tableRowColor), '@report_series.writeSeriesForTable: tableRowColor must be a string');
assert(isint(o.tableRowIndent), '@report_series.writeSeriesForTable: tableRowIndent must be an integer');
assert(islogical(o.tableShowMarkers), '@report_series.writeSeriesForTable: tableShowMarkers must be true or false');
assert(islogical(o.tableAlignRight), '@report_series.writeSeriesForTable: tableAlignRight must be true or false');
assert(isfloat(o.tableMarkerLimit), '@report_series.writeSeriesForTable: tableMarkerLimit must be a float');

%% Write Output
fprintf(fid, '%% Table Row (report_series)\n');
if ~isempty(o.tableRowColor)
    fprintf(fid, '\\rowcolor{%s}', o.tableRowColor);
end
if ~isempty(o.tableSubSectionHeader)
    fprintf(fid, '%s', o.tableSubSectionHeader);
    fprintf(fid, '\\\\%%\n');
    return;
end
if o.tableAlignRight
    fprintf(fid, '\\multicolumn{1}{r}{');
end
if o.tableRowIndent == 0
    fprintf(fid, '\\noindent');
else
    for i=1:o.tableRowIndent
        fprintf(fid,'\\indent');
    end
end
fprintf(fid, ' %s', o.data.tex{:});
if o.tableAlignRight
    fprintf(fid, '}');
end

printSeries(o, fid, o.data, dates{1}, precision);
if ~isempty(o.tableDataRhs)
    printSeries(o, fid, o.tableDataRhs, dates{2}, precision);
end

fprintf(fid, '\\\\%%\n');
end
