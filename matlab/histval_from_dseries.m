function [endo_histval, exo_histval] = histval_from_dseries(ds, initialperiod, DynareModel)

% Builds endo_histval and exo_hsitval from the content of a dseries object.
%
% INPUTS
% - ds                [dseries]    Dataset.
% - initialperiod     [dates]      Initial period of the simulation.
% - DynareModel       [struct]     Description of the model (M_ global structure).
%
% OUTPUTS
% - endo_histval      [double]     Vector of lagged values for the endogenous variables.
% - exo_histval       [double]     Matrix of lagged values for the exogenous variables.

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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

endo_histval = zeros(DynareModel.endo_nbr,DynareModel.maximum_endo_lag);

k = 1;
for i = 1:DynareModel.endo_nbr
    if i <= DynareModel.orig_endo_nbr
        if DynareModel.lead_lag_incidence(1,i) > 0
            if any(strcmp(deblank(DynareModel.endo_names(i,:)),ds.name))
                endo_histval(i,DynareModel.maximum_endo_lag) = ...
                    ds{deblank(DynareModel.endo_names(i,:))}(initialperiod-1).data;
            else
                error(sprintf('Can''t find %s in dseries', ...
                              deblank(DynareModel.endo_names(i,:))))
            end
        end
    else
        a = DynareModel.aux_vars(k);
        if a.type == 1
            if any(strcmp(deblank(DynareModel.endo_names(a.orig_index,:)), ds.name))
                endo_histval(i,DynareModel.maximum_endo_lag) = ...
                    ds{deblank(DynareModel.endo_names(a.orig_index,:))}(initialperiod-1+a.orig_lead_lag).data;
            else
                error(sprintf('Can''t find %s in dseries', ...
                              deblank(DynareModel.endo_names(a.orig_index,:))))
            end
        end
        k = k + 1;
    end
end

if nargout>1
    exo_histval = zeros(DynareModel.maximum_exo_lag, DynareModel.exo_nbr);
    exo_list = cellstr(DynareModel.exo_names);
    available_exo_variables = ismember(exo_list, ds.name);
    if any(~available_exo_variables)
        skipline()
        disp('Some exogenous variables are not available in the dseries object.')
        disp('Lagged values for these exogenous are zero.')
        skipline()
    end
    for t = 1:DynareModel.maximum_exo_lag
        for i=1:DynareModel.exo_nbr
            if available_exo_variables(i)
                exo_histval(DynareModel.maximum_exo_lag+1-t,i) = ds{exo_list{i}}(initialperiod-t).data;
            end
        end
    end
end
