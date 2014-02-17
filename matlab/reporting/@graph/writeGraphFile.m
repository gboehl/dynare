function o = writeGraphFile(o, pg, sec, row, col)
%function o = writeGraphFile(o, pg, sec, row, col)
% Write the tikz file that contains the graph
%
% INPUTS
%   o   [graph]   graph object
%   pg  [integer] this page number
%   sec [integer] this section number
%   row [integer] this row number
%   col [integer] this col number
%
% OUTPUTS
%   o   [graph] graph object
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

ne = length(o.series);
if ne < 1
    warning('@graph.writeGraphFile: no series to plot, returning');
    return;
end

if isempty(o.figname)
    o.figname = sprintf('%s/pg%d_sec%d_row%d_col%d.tex', o.figDirName, pg, sec, row, col);
end

[fid, msg] = fopen(o.figname, 'w');
if fid == -1
    error(['@graph.writeGraphFile: ' msg]);
end

fprintf(fid, '\\begin{tikzpicture}');

if isempty(o.xrange)
    dd = getMaxRange(o.series);
else
    dd = o.xrange;
end

fprintf(fid, '\\begin{axis}[%%\n');
% set tick labels
if isempty(o.xTickLabels)
    x = 1:1:dd.ndat;
    xTickLabels = strings(dd);
    fprintf(fid, 'xminorticks=true,\nyminorticks=true,\n');
else
    fprintf(fid,'minor xtick,\n');
    x = o.xTicks;
    xTickLabels = o.xTickLabels;
end
fprintf(fid, 'xticklabels={');
for i = 1:length(x)
    fprintf(fid,'%s,',lower(xTickLabels{i}));
end
fprintf(fid, '},\nxtick={');
for i = 1:length(x)
    fprintf(fid, '%d',x(i));
    if i ~= length(x)
        fprintf(fid,',');
    end
end
fprintf(fid, '},\nx tick label style={rotate=%f', o.xTickLabelRotation);
if o.xTickLabelRotation ~= 0
    fprintf(fid, ',anchor=%s', o.xTickLabelAnchor);
end
fprintf(fid, ['},\n',...
              'width=%fin,\n'...
              'height=%fin,\n'...
              'scale only axis,\n'...
              'axis lines=box,\n'...
              'unbounded coords=jump,\n'], o.width, o.height);

if ~isempty(o.title{1})
    fprintf(fid, 'title style={align=center');
    if ~isempty(o.titleFormat)
        fprintf(fid, ',font=%s', o.titleFormat);
    end
    fprintf(fid, '},\ntitle=');
    nt = length(o.title);
    for i=1:nt
        fprintf(fid, '%s', o.title{i});
        if i ~= nt
            fprintf(fid, '\\\\');
        end
    end
    fprintf(fid, ',\n');
end

if o.xAxisTight
    fprintf(fid, 'enlarge x limits=false,\n');
else
    fprintf(fid, 'enlarge x limits=true,\n');
end

if isempty(o.yrange)
    if o.yAxisTight
        fprintf(fid, 'enlarge y limits=false,\n');
    else
        fprintf(fid, 'enlarge y limits=true,\n');
    end
else
    fprintf(fid, 'ymin=%f,\nymax=%f,\n',o.yrange(1),o.yrange(2));
end

if o.showLegend
    fprintf(fid, 'legend style={');
    if ~o.showLegendBox
        fprintf(fid, 'draw=none,');
    end
    fprintf(fid, 'font=\\%s,', o.legendFontSize);
    if strcmp(o.legendOrientation, 'horizontal')
        fprintf(fid,'legend columns=-1,');
    end
    fprintf(fid, '},\nlegend pos=%s,\n', o.legendLocation);
end

if o.showGrid
    fprintf(fid, 'xmajorgrids=true,\nymajorgrids=true,\n');
end

if ~isempty(o.xlabel)
    fprintf(fid, 'xlabel=%s,\n', o.xlabel);
end

if ~isempty(o.ylabel)
    fprintf(fid, 'ylabel=%s,\n', o.ylabel);
end
fprintf(fid, ']\n');

if o.showZeroline
    fprintf(fid, '%%zeroline\n\\addplot[black,line width=.5,forget plot] coordinates {(1,0)(%d,0)};\n',dd.ndat);
end

for i=1:ne
    o.series{i}.writeSeriesForGraph(fid, dd);
    if o.showLegend
        fprintf(fid, '\\addlegendentry{%s}\n', o.series{i}.getTexName());
    end
end

if ~isempty(o.shade)
    xTickLabels = strings(dd);
    x1 = find(strcmpi(date2string(o.shade(1)), xTickLabels));
    x2 = find(strcmpi(date2string(o.shade(end)), xTickLabels));
    assert(~isempty(x1) && ~isempty(x2), ['@graph.writeGraphFile: either ' ...
                        date2string(o.shade(1)) ' or ' date2string(o.shade(end)) ' is not in the date ' ...
                        'range of data selected.']);
    fprintf(fid,['\\begin{pgfonlayer}{background}\n\\fill[%s!%f]\n(axis ' ...
                 'cs:%f,\\pgfkeysvalueof{/pgfplots/ymin})\nrectangle (axis ' ...
                 'cs:\\pgfkeysvalueof{/pgfplots/xmax},\\pgfkeysvalueof{/' ...
                 'pgfplots/ymax});\n\\end{pgfonlayer}\n'], ...
            o.shadeColor, o.shadeOpacity,x1);
end

fprintf(fid, '\\end{axis}\n\\end{tikzpicture}\n');
if fclose(fid) == -1
    error('@graph.writeGraphFile: closing %s\n', o.filename);
end
end