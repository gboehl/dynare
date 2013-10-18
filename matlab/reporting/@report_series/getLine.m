function h = getLine(o, xrange)
%function h = getLine(o, xrange)
% Create the report_series
%
% INPUTS
%   o       [report_series]    report_series object
%   xrange  [dynDates]  range of x values for line
%
% OUTPUTS
%   h       [handle]    handle to line
%
% SPECIAL REQUIREMENTS
%   none

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

%% Validate options provided by user
assert(~isempty(o.data) && isa(o.data, 'dynSeries'), ['@report_series.getLine: must ' ...
                    'provide data as a dynSeries']);

% Line
assert(ischar(o.graphLineColor), '@report_series.getLine: graphLineColor must be a string');
valid_line_style = {'none', '-', '--', ':', '-.'};
assert(any(strcmp(o.graphLineStyle, valid_line_style)), ...
       ['@report_series.getLine: graphLineStyle must be one of ' strjoin(valid_line_style, ' ')]);
assert(isfloat(o.graphLineWidth), ['@report_series.getLine: graphLineWidth must be a ' ...
                    'positive number']);

% GraphMarker
valid_graphMarker = {'+', 'o', '*', '.', 'x', 's', 'square', 'd', 'diamond', ...
                '^', 'v', '>', '<', 'p', 'pentagram', 'h', 'hexagram', ...
                'none'};
assert(isempty(o.graphMarker) || any(strcmp(o.graphMarker, valid_graphMarker)), ...
       ['@report_series.getLine: graphMarker must be one of ' strjoin(valid_graphMarker)]);

assert(ischar(o.graphMarkerEdgeColor), '@report_series.getLine: graphMarkerEdgeColor must be a string');
assert(ischar(o.graphMarkerFaceColor), '@report_series.getLine: graphMarkerFaceColor must be a string');
assert(isfloat(o.graphMarkerSize), ['@report_series.getLine: graphMarkerSize must be a ' ...
                    'positive number']);

% Marker & Line
assert(~(strcmp(o.graphLineStyle, 'none') && isempty(o.graphMarker)), ['@report_series.getLine: ' ...
                    'you must provide at least one of graphLineStyle and graphMarker']);

% Validate xrange
assert(isempty(xrange) || isa(xrange, 'dynDates'));

% Zero tolerance
assert(isfloat(o.zerotol), '@report_series.write: zerotol must be a float');

%%
if isempty(xrange) || xrange == o.data.time
    ds = o.data;
else
    ds = o.data(xrange);
end

% if graphing data that is within zerotol, set to zero, create report_series and
% get line:
thedata = ds.data;
stz = bsxfun(@and, ...
             bsxfun(@lt, thedata, o.zerotol), ...
             bsxfun(@gt, thedata, -o.zerotol));
if any(stz)
    thedata(stz) = 0;
end

opt = {'XData', 1:length(thedata)};
opt = {opt{:}, 'YData', thedata};

opt = {opt{:}, 'Color', o.graphLineColor};
opt = {opt{:}, 'LineStyle', o.graphLineStyle};
opt = {opt{:}, 'LineWidth', o.graphLineWidth};

if ~isempty(o.graphMarker)
    opt = {opt{:}, 'Marker', o.graphMarker};
    opt = {opt{:}, 'MarkerSize', o.graphMarkerSize};
    opt = {opt{:}, 'MarkerEdgeColor', o.graphMarkerEdgeColor};
    opt = {opt{:}, 'MarkerFaceColor', o.graphMarkerFaceColor};
end

h = line(opt{:});
end
