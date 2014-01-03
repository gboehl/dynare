function o = write(o, fid)
%function o = write(o, fid)
% Write a Report_Table object
%
% INPUTS
%   o           [report_table]    report_table object
%   fid         [integer]  file id
%
% OUTPUTS
%   o           [report_table]    report_table object
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

assert(fid ~= -1);
if ~o.seriesElements.numSeriesElements()
    warning('@report_table.write: no series to plot, returning');
    return;
end

%number of left-hand columns, 1 until we allow the user to group data,
% e.g.: GDP Europe
%         GDP France
%         GDP Germany
% this example would be two lh columns, with GDP Europe spanning both
nlhc = 1;

if isempty(o.range)
    dates = o.seriesElements.getMaxRange();
else
    dates = o.range{1};
end
ndates = dates.ndat;

fprintf(fid, '%% Report_Table Object\n');
fprintf(fid, '\\setlength{\\tabcolsep}{4pt}\n');
fprintf(fid, '\\begin{tabular}{@{}l');

for i=1:ndates
    if o.showVlines
        fprintf(fid, 'r|');
    else
        fprintf(fid, 'r');
        if o.vlineAfterEndOfPeriod
            if dates(i).time(2) == dates(i).freq
                fprintf(fid, '|');
            end
        end
        if ~isempty(o.vlineAfter)
            for j=1:length(o.vlineAfter)
                if dates(i) == o.vlineAfter{j}
                    if ~(o.vlineAfterEndOfPeriod && dates(i).time(2) == dates(i).freq)
                        fprintf(fid, '|');
                    end
                end
            end
        end
    end
end
datedata = dates.time;
years = unique(datedata(:, 1));
if length(o.range) > 1
    rhscols = strings(o.range{2});
    if o.range{2}.freq == 1
        rhscols = strrep(rhscols, 'Y', '');
    end
else
    rhscols = {};
end
for i=1:length(rhscols)
    fprintf(fid, 'r');
    if o.showVlines
        fprintf(fid, '|');
    end
end
nrhc = length(rhscols);
ncols = ndates+nlhc+nrhc;
fprintf(fid, '@{}}%%\n');
if ~isempty(o.title)
    fprintf(fid, '\\multicolumn{%d}{c}{\\%s %s}\\\\\n', ...
            ncols, o.titleSize, o.title);
end
fprintf(fid, '\\toprule%%\n');

% Column Headers
thdr = num2cell(years, size(years, 1));
switch dates.freq
    case 1
        for i=1:size(thdr, 1)
            fprintf(fid, ' & %d', thdr{i, 1});
        end
        for i=1:length(rhscols)
            fprintf(fid, ' & %s', rhscols{i});
        end
    case 4
        thdr{1, 2} = datedata(:, 2)';
        if size(thdr, 1) > 1
            for i=2:size(thdr, 1)
                split = find(thdr{i-1, 2} == 4, 1, 'first');
                assert(~isempty(split), '@report_table.write: Shouldn''t arrive here');
                thdr{i, 2} = thdr{i-1, 2}(split+1:end);
                thdr{i-1, 2} = thdr{i-1, 2}(1:split);
            end
        end
        for i=1:size(thdr, 1)
            fprintf(fid, ' & \\multicolumn{%d}{c}{%d}', size(thdr{i,2}, 2), thdr{i,1});
        end
        for i=1:length(rhscols)
            fprintf(fid, ' & %s', rhscols{i});
        end
        fprintf(fid, '\\\\\\cline{%d-%d}%%\n', nlhc+1, ncols);
        for i=1:size(thdr, 1)
            quarters = thdr{i, 2};
            for j=1:size(quarters, 2)
                fprintf(fid, ' & \\multicolumn{1}{c}{Q%d}', quarters(j));
            end
        end
    case 12
        thdr{1, 2} = datedata(:, 2)';
        if size(thdr, 1) > 1
            for i=2:size(thdr, 1)
                split = find(thdr{i-1, 2} == 12, 1, 'first');
                assert(~isempty(split), '@report_table.write 2: Shouldn''t arrive here');
                thdr{i, 2} = thdr{i-1, 2}(split+1:end);
                thdr{i-1, 2} = thdr{i-1, 2}(1:split);
            end
        end
        for i=1:size(thdr, 1)
            fprintf(fid, ' & \\multicolumn{%d}{c}{%d}', size(thdr{i,2}, 2), thdr{i,1});
        end
        for i=1:length(rhscols)
            fprintf(fid, ' & %s', rhscols{i});
        end
        fprintf(fid, '\\\\\\cline{%d-%d}%%\n', nlhc+1, ncols);
        for i=1:size(thdr, 1)
            months = thdr{i, 2};
            for j=1:size(months, 2)
                fprintf(fid, ' & \\multicolumn{1}{c}{M%d}', months(j));
            end
        end
    otherwise
        error('@report_table.write: invalid dseries frequency');
end
fprintf(fid, '\\\\[-2pt]%%\n');
fprintf(fid, '\\hline%%\n');
fprintf(fid, '%%\n');

% Write Report_Table Data
ne = o.seriesElements.numSeriesElements();
for i=1:ne
    o.seriesElements(i).write(fid, o.range, o.precision);
    if o.showHlines
        fprintf(fid, '\\hline\n');
    end
end

fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular} \\par \\medskip\n\n');
fprintf(fid, '%% End Report_Table Object\n');
end
