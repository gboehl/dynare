function set_historical_values(ds, initialperiod)

% Builds endo_histval and exo_hsitval from the content of a dseries object.
%
% INPUTS
% - ds                [dseries]    Dataset.
% - initialperiod     [dates]      Initial period of the simulation.
%
% OUTPUTS
% - none

% Copyright (C) 2018 Dynare Team
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

global M_

if ischar(ds)
    ds = evalin('caller', ds);
end

if ischar(initialperiod)
    initialperiod = eval(initialperiod);
end

% Initialize endo_histval.
M_.endo_histval = zeros(M_.endo_nbr, M_.orig_maximum_endo_lag);

% Fill endo_histval.
k = 1;
for i = 1:M_.endo_nbr
    if i <= M_.orig_endo_nbr
        if M_.lead_lag_incidence(1,i) > 0
            if any(strcmp(M_.endo_names{i}, ds.name))
                M_.endo_histval(i, M_.maximum_endo_lag) = ...
                    ds{M_.endo_names{i}}(initialperiod).data;
            else
                error(sprintf('Can''t find %s in dseries', M_.endo_names{i}))
            end
        end
    else
        a = M_.aux_vars(k);
        if a.type == 1
            if any(strcmp(M_.endo_names{a.orig_index}, ds.name))
                M_.endo_histval(i,M_.maximum_endo_lag) = ...
                    ds{M_.endo_names{a.orig_index}}(initialperiod+a.orig_lead_lag).data;
            else
                error(sprintf('Can''t find %s in dseries', M_.endo_names{a.orig_index}))
            end
        end
        k = k + 1;
    end
end

% If model has lags on exogenous variables, initialize and fill exo_histval
if M_.maximum_exo_lag
    M_.exo_histval = zeros(M_.maximum_exo_lag, M_.exo_nbr);
    exo_list = M_.exo_names;
    available_exo_variables = ismember(exo_list, ds.name);
    if any(~available_exo_variables)
        skipline()
        disp('Some exogenous variables are not available in the dseries object.')
        disp('Default value for lagged exogenous variables is zero.')
        skipline()
    end
    for t = 1:M_.maximum_exo_lag
        for i=1:M_.exo_nbr
            if available_exo_variables(i)
                M_.exo_histval(M_.maximum_exo_lag+1-t,i) = ds{exo_list{i}}(initialperiod-t).data;
            end
        end
    end
end