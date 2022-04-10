function [endogenousvariables, info] = sim1_purely_backward(endogenousvariables, exogenousvariables, steadystate, M, options)

% Performs deterministic simulation of a purely backward model

% Copyright © 2012-2022 Dynare Team
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

ny0 = nnz(M.lead_lag_incidence(2,:));    % Number of variables at current period
iyb = M.lead_lag_incidence(1,:)>0;       % Logical vector (for lagged variables)

if ny0 ~= M.endo_nbr
    error('All endogenous variables must appear at the current period!')
end

if ismember(options.solve_algo, [12,14]) && ~M.possible_to_use_solve_algo_12_14
    error(M.message_solve_algo_12_14)
end

dynamicmodel = str2func(sprintf('%s.%s', M.fname, 'dynamic'));
dynamicmodel_s = str2func('dynamic_backward_model_for_simulation');

info.status = true;

for it = M.maximum_lag + (1:options.periods)
    y = endogenousvariables(:,it-1);        % Values at previous period, also used as guess value for current period
    ylag = y(iyb);
    if ismember(options.solve_algo, [12,14])
        [tmp, check, ~, ~, errorcode] = dynare_solve(dynamicmodel_s, y, ...
                                                     options.simul.maxit, options.dynatol.f, options.dynatol.x, ...
                                                     options, M.isloggedlhs, M.isauxdiffloggedrhs, M.endo_names, M.lhs, ...
                                                     dynamicmodel, ylag, exogenousvariables, M.params, steadystate, it);
    else
        [tmp, check, ~, ~, errorcode] = dynare_solve(dynamicmodel_s, y, ...
                                                     options.simul.maxit, options.dynatol.f, options.dynatol.x, ...
                                                     options, dynamicmodel, ylag, exogenousvariables, M.params, steadystate, it);
    end
    if check
        info.status = false;
        dprintf('sim1_purely_backward: Nonlinear solver routine failed with errorcode=%i. in period %i', errorcode, it)
        break
    end
    endogenousvariables(:,it) = tmp;
end