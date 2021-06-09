function [endogenousvariables, info] = sim1_purely_static(endogenousvariables, exogenousvariables, steadystate, M, options)

% Performs deterministic simulation of a purely static model

% Copyright © 2021 Dynare Team
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

if nnz(M.lead_lag_incidence(1,:)) ~= M.endo_nbr
    error('All endogenous variables must appear at the current period!')
end

if ismember(options.solve_algo, [12,14]) && ~M.possible_to_use_solve_algo_12_14
    error(M.message_solve_algo_12_14)
end

dynamicmodel = str2func([M.fname,'.dynamic']);
dynamicmodel_s = str2func('dynamic_static_model_for_simulation');

info.status = true;

y = endogenousvariables(:,1);

for it = 1:options.periods
    if ismember(options.solve_algo, [12,14])
        [tmp, check] = dynare_solve(dynamicmodel_s, y, options, M.isloggedlhs, M.isauxdiffloggedrhs, M.endo_names, M.lhs, ...
                                    dynamicmodel, exogenousvariables, M.params, steadystate, it);
    else
        [tmp, check] = dynare_solve(dynamicmodel_s, y, options, ...
                                    dynamicmodel, exogenousvariables, M.params, steadystate, it);
    end
    if check
        info.status = false;
    end
    endogenousvariables(:,it) = tmp;
    y = endogenousvariables(:,it);
end