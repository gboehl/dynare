function [endogenousvariables, info] = sim1_purely_forward(endogenousvariables, exogenousvariables, steadystate, M, options)
% Performs deterministic simulation of a purely forward model

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

ny0 = nnz(M.lead_lag_incidence(1,:));    % Number of variables at current period
iyf = find(M.lead_lag_incidence(2,:)>0); % Indices of variables at next period

if ny0 ~= M.endo_nbr
    error('All endogenous variables must appear at the current period!')
end

dynamicmodel = str2func(sprintf('%s.%s', M.fname, 'dynamic'));
dynamicmodel_s = str2func('dynamic_forward_model_for_simulation');

info.status = 1;

for it = options.periods:-1:1
    yf = endogenousvariables(:,it+1); % Values at next period, also used as guess value for current period
    if ismember(options.solve_algo, [12,14])
        [tmp, check, ~, ~, errorcode] = dynare_solve(dynamicmodel_s, yf, options, M.isloggedlhs, M.isauxdiffloggedrhs, M.endo_names, M.lhs, ...
                                                     dynamicmodel, yf(iyf), exogenousvariables, M.params, steadystate, it);
    else
        [tmp, check, ~, ~, errorcode] = dynare_solve(dynamicmodel_s, yf, options, ...
                                                     dynamicmodel, yf(iyf), exogenousvariables, M.params, steadystate, it);
    end
    if check
        info.status = 0;
        if option.debug
            dprintf('sim1_purely_forward: Nonlinear solver routine failed with errorcode=%i in period %i.', errorcode, it)
        end
    end
    endogenousvariables(:,it) = tmp(1:M.endo_nbr);
end
