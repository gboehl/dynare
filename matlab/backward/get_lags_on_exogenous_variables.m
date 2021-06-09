function l = get_lags_on_exogenous_variables(DynareModel)

% Returns a vector with the max lag for each exogenous variable.

% Copyright (C) 2017 Dynare Team
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

l = zeros(DynareModel.exo_nbr, 1);

if ~isempty(DynareModel.aux_vars)
    aux_var_for_lagged_exogenous = find([DynareModel.aux_vars(:).type]==3);
    for i=1:length(aux_var_for_lagged_exogenous)
        l(DynareModel.aux_vars(aux_var_for_lagged_exogenous(i)).orig_index) = ...
            DynareModel.aux_vars(aux_var_for_lagged_exogenous(i)).orig_lead_lag-1;
    end
end