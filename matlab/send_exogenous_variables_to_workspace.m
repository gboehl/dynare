function send_exogenous_variables_to_workspace()
% send_exogenous_variables_to_workspace()
% Saves all the endogenous variables in matlab's workspace.

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
global M_ oo_

for idx = 1:M_.exo_nbr
    assignin('base', M_.exo_names{idx}, oo_.exo_simul(:,idx)')
end