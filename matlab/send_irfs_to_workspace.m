function send_irfs_to_workspace()
% Saves all the IRFs in MATLAB's workspace.

% Copyright © 2023 Dynare Team
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
global oo_

if isfield(oo_,'irfs')
    irf_fields=fieldnames(oo_.irfs);
    for irf_iter = 1:size(irf_fields,1)
    assignin('base',irf_fields{irf_iter},oo_.irfs.(irf_fields{irf_iter})');
    end
end
