function str = subst_auxvar(var_index, aux_lead_lag, M_)
% Given the index of an endogenous (possibly an auxiliary var), and a
% lead/lag, creates a string of the form "x(lag)".
% In the case of auxiliary vars for lags, replace by the original variable
% name, and compute the lead/lag accordingly.
% INPUTS
% - aux_index           [double]    index of auxiliary variable
% - aux_lead_lag        [double]    index for lead or lag
% - M_                  [struct]    model structure
%
% OUTPUTS
% - str                 [string]    name of auxiliary

% Copyright (C) 2001-2021 Dynare Team
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

if var_index <= M_.orig_endo_nbr
    str = sprintf('%s(%d)', M_.endo_names{var_index}, aux_lead_lag);
    return
end
aux_index=find([M_.aux_vars(:).endo_index]==var_index);
if ~isempty(aux_index)
    switch M_.aux_vars(aux_index).type
        case 0
            % endo leads >= 2
            str = sprintf('%s(%d)', M_.endo_names{var_index}, aux_lead_lag);
            return
        case 1
            % endo lags >= 2
            orig_name = M_.endo_names{M_.aux_vars(aux_index).orig_index};
        case 2
            % exo leads >= 1
            orig_name = M_.exo_names{M_.aux_vars(aux_index).orig_index};
        case 3
            % exo lags >= 1
            orig_name = M_.exo_names{M_.aux_vars(aux_index).orig_index};
        case 4
            % Expectation operator
            str = sprintf('EXPECTATION(%d)(...)', aux_lead_lag);
            return
        case 6
            % Ramsey's multipliers
            if ~isempty(aux_lead_lag)
                str = sprintf('%s(%d)', M_.endo_names{M_.aux_vars(aux_index).endo_index}, aux_lead_lag);
            else
                str = sprintf('%s', M_.endo_names{M_.aux_vars(aux_index).endo_index});
            end
            return
        case 7
            % currently unused
            error('This type of auxiliary variable should not occur, please contact the developers.')
        case 8
            % Diff operator
            str = sprintf('diff(%s)', M_.endo_names{M_.aux_vars(aux_index).orig_index});
            return
        case 9
            % Lagged diff
            lags = 0;
            j = aux_index;
            while M_.aux_vars(j).type==9
                j = M_.aux_vars(j).orig_index;
                lags = lags+1;
            end
            str = sprintf('diff(%s(-%u))', M_.endo_names{M_.aux_vars(j).orig_index}, lags);
            return
        case 10
            % Variable created when diff was taken of unary operator (log, exp)
            error('This type of auxiliary variable should not occur, please contact the developers.')
        case 11
            % Leaded diff
            leads = 0;
            j = aux_index;
            while M_.aux_vars(j).type==11
                j = M_.aux_vars(j).orig_index;
                lags = lags+1;
            end
            str = sprintf('diff(%s(%u))', M_.endo_names{M_.aux_vars(j).orig_index}, leads);
            return
        otherwise
            error('Invalid auxiliary type: %s', M_.endo_names{var_index})
    end
    str = sprintf('%s(%d)', orig_name, M_.aux_vars(aux_index).orig_lead_lag+aux_lead_lag);
    return
else
    error('Could not find aux var: %s', M_.endo_names{var_index})
end
