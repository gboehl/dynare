function l = get_lags_on_endogenous_variables(DynareModel)

% Returns a vector with the max lag for each endogenous variable.

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

l = zeros(DynareModel.orig_endo_nbr, 1);
l(find(DynareModel.lead_lag_incidence(1,1:DynareModel.orig_endo_nbr))) = -1;

if ~isempty(DynareModel.aux_vars)
    aux_var_for_lagged_endogenous = find([DynareModel.aux_vars(:).type]==1);
    for i=1:length(aux_var_for_lagged_endogenous)
        l(DynareModel.aux_vars(aux_var_for_lagged_endogenous(i)).orig_index) = ...
            DynareModel.aux_vars(aux_var_for_lagged_endogenous(i)).orig_lead_lag;
    end
end