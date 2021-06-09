function DynareModel = set_exogenous_variables_for_simulation(DynareModel)

% Appends the list of observed exogenous variables in Dynare's model structure (if any).
%
% INPUTS
% - DynareModel   [struct]    Dynare's model global structure, M_.
%
% OUTPUTS
% - DynareModel   [struct]    Dynare's model global structure, M_.
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2019 Dynare Team
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

if isfield(DynareModel, 'exo_partitions')
    if isfield(DynareModel.exo_partitions, 'used')
        DynareModel.simulation_exo_names = DynareModel.exo_names(~strcmpi('estimationonly', DynareModel.exo_partitions.used));
    end
end