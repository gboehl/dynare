function o = createGraph(o)
%function o = createGraph(o)
% Create the graph
%
% INPUTS
%   o   [graph] graph object
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

if ~isempty(o.figname)
    warning('@graph.createGraph: will overwrite %s with new graph\n', ...
            o.figname);
end

if ~o.seriesElements.numSeriesElements()
    warning('@graph.crepateGraph: no series to plot, returning');
    return;
end

if isempty(o.graphSize)
    h = figure('visible','off');
else
    h = figure('visible','off','position',[1, 1, o.graphSize(1), o.graphSize(2)]);
end
hold on;
box on;
if o.showGrid
    grid on;
    set(gca, 'GridLineStyle', '--');
end

if isempty(o.xrange)
    dd = o.seriesElements.getMaxRange();
else
    dd = o.xrange;
end

ne = o.seriesElements.numSeriesElements();
line_handles = zeros(ne, 1);
for i=1:ne
    line_handles(i) = o.seriesElements(i).getLine(dd);
end

x = 1:1:dd.ndat;
xlim([1 dd.ndat]);
xlabels = strings(dd);

if ~isempty(o.yrange)
    ylim(o.yrange);
end

if o.showZeroline
    a = ylim;
    if 0 > a(1) && 0 < a(2)
        lh = line(xlim, [0 0], 'color', 'k', 'LineWidth', 0.25);
        children =get(gca(), 'children');
        children = [children(2:end); lh];
        set(gca(), 'children', children);
    end
end

if ~isempty(o.shade)
    x1 = find(strcmpi(date2string(o.shade(1)), xlabels));
    x2 = find(strcmpi(date2string(o.shade(end)), xlabels));
    assert(~isempty(x1) && ~isempty(x2), ['@graph.createGraph: either ' ...
                        date2string(o.shade(1)) ' or ' date2string(o.shade(end)) ' is not in the date ' ...
                        'range of data selected.']);
    yrange = get(gca, 'YLim');

    % From ShadePlotForEmpahsis (Matlab Exchange)
    % use patch bc area doesn't work with matlab2tikz
    sh = patch([repmat(x1, 1, 2) repmat(x2, 1, 2)], ...
               [yrange fliplr(yrange)], o.shadeColor, ...
               'facealpha', o.shadeOpacity);
    children = get(gca, 'children');
    children = [children(2:end); sh];
    set(gca, 'children', children);
end

if isempty(o.xTickLabels)
    xticks = get(gca, 'XTick');
    xTickLabels = cell(1, length(xticks));
    for i=1:length(xticks)
        if xticks(i) >= x(1) && ...
                xticks(i) <= x(end)
            xTickLabels{i} = xlabels{xticks(i)};
        else
            xTickLabels{i} = '';
        end
    end
else
    set(gca, 'XTick', o.xTicks);
    xTickLabels = o.xTickLabels;
end
set(gca, 'XTickLabel', xTickLabels);

if o.showLegend
    lh = legend(line_handles, o.seriesElements.getTexNames(), ...
                'orientation', o.legendOrientation, ...
                'location', o.legendLocation);
    set(lh, 'FontSize', o.legendFontSize);
    set(lh, 'interpreter', 'latex');
    if ~o.showLegendBox
        legend('boxoff');
    end
end

if ~isempty(o.xlabel)
    xlabel(['$\textbf{\footnotesize ' o.xlabel '}$'], 'Interpreter', 'LaTex');
end

if ~isempty(o.ylabel)
    ylabel(['$\textbf{\footnotesize ' o.ylabel '}$'], 'Interpreter', 'LaTex');
end
drawnow;

if isempty(o.figname)
    [junk, tn] = fileparts(tempname);
    if strcmp(computer, 'PCWIN') || strcmp(computer, 'PCWIN64')
        tn = strrep(tn, '_', '-');
    end
    o.figname = [o.figDirName '/' tn '.tex'];
end

if isoctave && isempty(regexpi(computer, '.*apple.*', 'once'))
    print(o.figname, '-dtikz');
else
    matlab2tikz('filename', o.figname, ...
                'showInfo', false, ...
                'showWarnings', false, ...
                'checkForUpdates', false);
end

grid off;
box off;
hold off;
close(h);
clear h;
fprintf(1, '.');
end
