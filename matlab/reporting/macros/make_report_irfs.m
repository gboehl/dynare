function make_report_irfs(M_, oo_, ticks_every, showOutput)
% Builds canned IRF report
%
% INPUTS
%   M_            [struct]
%   oo_           [struct]
%   ticks_every   [int]      number of spaces between ticks. Default 5.
%   showOutput    [bool]     showOutput the report. Default true
%
% OUTPUTS
%   None
%
% SPECIAL REQUIREMENTS
%   None

% Copyright © 2015-2023 Dynare Team
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

if ~isfield(oo_, 'irfs')
    disp('make_report_irfs: oo_.irfs does not exist');
    return
end
fields = fieldnames(oo_.irfs);
if isempty(fields)
    disp('make_report_irfs: oo_.irfs is empty');
    return
end
if ~isfield(M_, 'exo_names')
    disp('make_report_irfs: M_.exo_names does not exist');
    return
end
if ~isfield(M_, 'endo_names')
    disp('make_report_irfs: M_.endo_names does not exist');
    return
end
if ~isfield(M_, 'fname')
    disp('make_report_irfs: M_.fname does not exist');
    return
end

if nargin < 3
    ticks_every = 5;
else
    assert(isint(ticks_every));
end
if nargin < 4
    showOutput = true;
else
    assert(islogical(showOutput));
end
n6 = 1;
justAddedPage = 0;
calcxticks = false;
r = report('filename', [M_.fname '_canned_irf_report.tex'], 'showOutput', showOutput);
for i = 1:length(M_.exo_names)
    newexo = 1;
    for j = 1:length(M_.endo_names)
        idx = ismember(fields, [M_.endo_names{j} '_' M_.exo_names{i}]);
        if (mod(n6 - 1, 6) == 0 && ~justAddedPage) || ...
                (newexo && any(idx))
            r = r.addPage('title', {'Canned Irf Report'; ['shock ' ...
                                strrep(M_.exo_names{i},'_','\_')]});
            r = r.addSection('cols', 2);
            n6 = 1;
            justAddedPage = 1;
            newexo = 0;
        end
        if any(idx)
            if ~calcxticks
                data = dseries(oo_.irfs.(fields{idx})');
                xTicks = 1:ticks_every:floor(data.nobs/ticks_every)*ticks_every+1;
                xTickLabels = regexp(num2str(xTicks-1), '(?:\s)+', 'split');
                calcxticks = true;
            end
            r = r.addGraph('data', dseries(oo_.irfs.(fields{idx})'), ...
                           'title', M_.endo_names{j}, '_', '\_'), ...
                           'titleFormat', '\Huge', ...
                           'showGrid', false, ...
                           'yTickLabelZeroFill', false, ...
                           'yTickLabelPrecision', 1, ...
                           'showZeroLine', true, ...
                           'zeroLineColor', 'red', ...
                           'xTicks', xTicks, ...
                           'xTickLabels', xTickLabels);
            n6 = n6 + 1;
            justAddedPage = 0;
        end
    end
end
r.write();
r.compile();
end
