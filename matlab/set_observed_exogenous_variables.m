function M_ = set_observed_exogenous_variables(M_)
% M_ = set_observed_exogenous_variables(M_)
% Appends the list of observed exogenous variables in Dynare's model structure (if any).
%
% INPUTS
% - M_   [struct]    Dynare's model global structure.
%
% OUTPUTS
% - M_   [struct]    Dynare's model global structure.
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2019-2023 Dynare Team
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

if isfield(M_, 'exo_partitions')
    if isfield(M_.exo_partitions, 'status')
        M_.observed_exo_names = M_.exo_names(strcmpi('observed', M_.exo_partitions.status));
    end
end