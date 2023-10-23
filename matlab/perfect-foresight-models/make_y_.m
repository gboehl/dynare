function oo_=make_y_(M_, options_, oo_)

% forms oo_.endo_simul as guess values for deterministic simulations
%
% INPUTS
% - M_          [struct]   Dynare model structure
% - options_    [struct]   Dynare options structure
% - oo_         [struct]   Dynare results structure
%
% OUTPUTS
% - oo_         [struct]   Updated dynare results structure

% Copyright © 1996-2023 Dynare Team
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

if isempty(oo_.steady_state)
    oo_.steady_state = zeros(M_.endo_nbr,1);
end

if isempty(oo_.initval_series)
    if isempty(M_.endo_histval)
        if isempty(oo_.initial_steady_state)
            oo_.endo_simul = repmat(oo_.steady_state, 1, M_.maximum_lag+options_.periods+M_.maximum_lead);
        else
            oo_.endo_simul = [repmat(oo_.initial_steady_state, 1, M_.maximum_lag) repmat(oo_.steady_state, 1,options_.periods+M_.maximum_lead)];
        end
    else
        if ~isempty(oo_.initial_steady_state)
            error('histval and endval cannot be used simultaneously')
        end
        % the first NaNs take care of the case where there are lags > 1 on
        % exogenous variables
        oo_.endo_simul = [M_.endo_histval ...
                          repmat(oo_.steady_state, 1, options_.periods+M_.maximum_lead)];
    end
else
    y = oo_.initval_series{M_.endo_names{:}}.data;
    oo_.endo_simul = y(M_.orig_maximum_lag - M_.maximum_lag + 1:M_.orig_maximum_lag + options_.periods + ...
                       M_.maximum_lead, :)';
    if ~isempty(M_.endo_histval)
        if ~isempty(oo_.initial_steady_state)
            error('histval and endval cannot be used simultaneously')
        end
        oo_.endo_simul(:,1:M_.maximum_lag) ...
            = M_.endo_histval(:, 1:M_.maximum_lag);
    end
end
