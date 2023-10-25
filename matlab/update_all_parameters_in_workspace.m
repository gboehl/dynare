function update_all_parameters_in_workspace(M_)
% update_all_parameters_in_workspace(M_)
% Updates all parameter values in Matlab/Octave base workspace.

% Copyright © 2018-2023 Dynare Team
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

for i=1:length(M_.params)
    assignin('base', M_.param_names{i}, M_.params(i));
end