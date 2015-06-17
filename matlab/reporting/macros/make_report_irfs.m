function make_report_irfs(oo)
% Builds posterior IRFs after the MH algorithm.
%
% INPUTS
%   oo            [struct]
%
% OUTPUTS
%   None
%
% SPECIAL REQUIREMENTS
%   None

% Copyright (C) 2015 Dynare Team
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

  if ~isfield(oo, 'irfs')
    disp('make_report_irfs: oo_.irfs does not exist');
    return
  end
  fields = fieldnames(oo.irfs);
  if isempty(fields)
    disp('make_report_irfs: oo_.irfs is empty');
    return
  end

  r = report();
  for i = 1:length(fields)
    if mod(i-1, 6) == 0
      r = r.addPage('title', {'Canned Irf Report'});
      r = r.addSection('cols', 2);
    end
    r = r.addGraph('data', dseries(oo.irfs.(fields{i})'), ...
                   'title', strrep(fields{i}, '_', '\_'), ...
                   'showGrid', false, ...
                   'yTickLabelZeroFill', false, ...
                   'showZeroLine', true, ...
                   'zeroLineColor', 'red');
  end
  r.write();
  r.compile();
end
