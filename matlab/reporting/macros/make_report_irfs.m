function make_report_irfs(M, oo)
% Builds canned IRF report
%
% INPUTS
%   M             [struct]
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
  if ~isfield(M, 'exo_names')
      disp('make_report_irfs: M_.exo_names does not exist');
      return
  end
  if ~isfield(M, 'endo_names')
      disp('make_report_irfs: M_.endo_names does not exist');
      return
  end
  
  n6 = 1;
  justAddedPage = 0;
  r = report();
  for i = 1:length(M.exo_names)
      for j = 1:length(M.endo_names)
          if mod(n6 - 1, 6) == 0 && ~justAddedPage
              r = r.addPage('title', {'Canned Irf Report'; ['shock ' ...
                  strrep(strtrim(M.exo_names(i,:)),'_','\_')]});
              r = r.addSection('cols', 2);
              n6 = 1;
              justAddedPage = 1;
          end
          idx = ismember(fields,[strtrim(M.endo_names(j,:)) '_' ...
              strtrim(M.exo_names(i,:))]);
          if any(idx)
              r = r.addGraph('data', dseries(oo.irfs.(fields{idx})'), ...
                  'title', strrep(M.endo_names(j,:), '_', '\_'), ...
                  'titleFormat', '\Huge', ...
                  'showGrid', false, ...
                  'yTickLabelZeroFill', false, ...
                  'showZeroLine', true, ...
                  'zeroLineColor', 'red');
              n6 = n6 + 1;
              justAddedPage = 0;
          end
      end
  end
  r.write();
  r.compile();
end
