function o = writeSeriesForGraph(o, fid, xrange)
%function o = writeSeriesForGraph(o, fid, xrange)
% Print a TikZ line
%
% INPUTS
%   o       [report_series]    series object
%   xrange  [dates]            range of x values for line
%
% OUTPUTS
%   NONE
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2014 Dynare Team
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
assert(~isempty(o.data) && isa(o.data, 'dseries'), ['@report_series.writeSeriesForGraph: must ' ...
                    'provide data as a dseries']);

% Line
valid_graphLineColor = {'red', 'green', 'blue', 'cyan ', 'magenta', 'yellow', ...
                    'black', 'gray', 'darkgray', 'lightgray', 'brown', ...
                    'lime', 'olive', 'orange', 'pink', 'purple', 'teal', 'violet', 'white'};
assert(any(strcmp(o.graphLineColor, valid_graphLineColor)), ...
       ['@report_series.writeSeriesForGraph: graphLineColor must be one of ' strjoin(valid_graphLineColor)]);
assert(ischar(o.graphLineStyle), '@report_series.writeSeriesForGraph: graphLineStyle must be a string');
assert(isfloat(o.graphLineWidth) && o.graphLineWidth > 0, ...
                    '@report_series.writeSeriesForGraph: graphLineWidth must be a positive number');

% GraphMarker
valid_graphMarker = {'x', '+', '-', '|', 'o', 'asterisk', 'star', '10-pointed star', 'oplus', ...
                    'oplus*', 'otimes', 'otimes*', 'square', 'square*', 'triangle', 'triangle*', 'diamond', ...
                    'diamond*', 'halfdiamond*', 'halfsquare*', 'halfsquare right*', ...
                    'halfsquare left*','Mercedes star','Mercedes star flipped','halfcircle',...
                    'halfcircle*','pentagon','pentagon star'};
assert(isempty(o.graphMarker) || any(strcmp(o.graphMarker, valid_graphMarker)), ...
       ['@report_series.writeSeriesForGraph: graphMarker must be one of ' strjoin(valid_graphMarker)]);

assert(ischar(o.graphMarkerEdgeColor), '@report_series.writeSeriesForGraph: graphMarkerEdgeColor must be a string');
assert(ischar(o.graphMarkerFaceColor), '@report_series.writeSeriesForGraph: graphMarkerFaceColor must be a string');
assert(isfloat(o.graphMarkerSize) && o.graphMarkerSize > 0, ...
                    '@report_series.writeSeriesForGraph: graphMarkerSize must be a positive number');

% Marker & Line
assert(~(strcmp(o.graphLineStyle, 'none') && isempty(o.graphMarker)), ['@report_series.writeSeriesForGraph: ' ...
                    'you must provide at least one of graphLineStyle and graphMarker']);

% Validate xrange


% Zero tolerance
assert(isfloat(o.zeroTol), '@report_series.write: zeroTol must be a float');

%%
if isempty(xrange) || all(xrange == o.data.dates)
    ds = o.data;
else
    ds = o.data(xrange);
end

% if graphing data that is within zeroTol, set to zero, create report_series and
% get line:
thedata = ds.data;
stz = bsxfun(@and, ...
             bsxfun(@lt, thedata, o.zeroTol), ...
             bsxfun(@gt, thedata, -o.zeroTol));
if any(stz)
    thedata(stz) = 0;
end

fprintf(fid, '%%series %s\n\\addplot[color=%s,%s,line width=%fpt,line join=round',...
        o.data.name{:}, o.graphLineColor, o.graphLineStyle, o.graphLineWidth);
if ~isempty(o.graphMarker)
    if isempty(o.graphMarkerEdgeColor)
        o.graphMarkerEdgeColor = o.graphLineColor;
    end
    if isempty(o.graphMarkerFaceColor)
        o.graphMarkerFaceColor = o.graphLineColor;
    end
    fprintf(fid, ',mark=%s,mark size=%f,every mark/.append style={draw=%s,fill=%s}',...
            o.graphMarker,o.graphMarkerSize,o.graphMarkerEdgeColor,o.graphMarkerFaceColor);
end
fprintf(fid,']\ntable[row sep=crcr]{\nx y\\\\\n');
for i=1:ds.dates.ndat
    fprintf(fid, '%d %f\\\\\n', i, thedata(i));
end
fprintf(fid,'};\n');
end
