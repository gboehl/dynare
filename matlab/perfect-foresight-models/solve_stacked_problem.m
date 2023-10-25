function [endogenousvariables, success, maxerror] = solve_stacked_problem(endogenousvariables, exogenousvariables, steadystate, M_, options_)
% [endogenousvariables, success, maxerror] = solve_stacked_problem(endogenousvariables, exogenousvariables, steadystate, M_, options_)
% Solves the perfect foresight model using dynare_solve
%
% INPUTS
% - endogenousvariables [double] N*T array, paths for the endogenous variables (initial guess).
% - exogenousvariables  [double] T*M array, paths for the exogenous variables.
% - steadystate         [double] N*1 array, steady state for the endogenous variables.
% - M_                   [struct] contains a description of the model.
% - options_             [struct] contains various options.
%
% OUTPUTS
% - endogenousvariables [double] N*T array, paths for the endogenous variables (solution of the perfect foresight model).
% - success             [logical] Whether a solution was found
% - maxerror            [double] 1-norm of the residual

% Copyright © 2015-2023 Dynare Team
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

[options_, y0, yT, z, i_cols, i_cols_J1, i_cols_T, i_cols_j, i_cols_1, i_cols_0, i_cols_J0, dynamicmodel] = ...
    initialize_stacked_problem(endogenousvariables, options_, M_, steadystate);

if (options_.solve_algo == 10 || options_.solve_algo == 11)% mixed complementarity problem
    [lb,ub,eq_index] = get_complementarity_conditions(M_,options_.ramsey_policy);
    if options_.linear_approximation
        lb = lb - steadystate_y;
        ub = ub - steadystate_y;
    end
    if options_.solve_algo == 10
        options_.lmmcp.lb = repmat(lb,options_.periods,1);
        options_.lmmcp.ub = repmat(ub,options_.periods,1);
    elseif options_.solve_algo == 11
        options_.mcppath.lb = repmat(lb,options_.periods,1);
        options_.mcppath.ub = repmat(ub,options_.periods,1);
    end
    [y, check, res, ~, errorcode] = dynare_solve(@perfect_foresight_mcp_problem, z(:), ...
                                               options_.simul.maxit, options_.dynatol.f, options_.dynatol.x, ...
                                               options_, ...
                                               dynamicmodel, y0, yT, ...
                                               exogenousvariables, M_.params, steadystate, ...
                                               M_.maximum_lag, options_.periods, M_.endo_nbr, i_cols, ...
                                               i_cols_J1, i_cols_1, i_cols_T, i_cols_j, i_cols_0, i_cols_J0, ...
                                               eq_index);
    eq_to_ignore=find(isfinite(lb) | isfinite(ub));

else
    [y, check, res, ~, errorcode] = dynare_solve(@perfect_foresight_problem, z(:), ...
                                               options_.simul.maxit, options_.dynatol.f, options_.dynatol.x, ...
                                               options_, y0, yT, exogenousvariables, M_.params, steadystate, options_.periods, M_, options_);
end

if all(imag(y)<.1*options_.dynatol.x)
    if ~isreal(y)
        y = real(y);
    end
else
    check = 1;
end

endogenousvariables(:, M_.maximum_lag+(1:options_.periods)) = reshape(y, M_.endo_nbr, options_.periods);
residuals=zeros(size(endogenousvariables));
residuals(:, M_.maximum_lag+(1:options_.periods)) = reshape(res, M_.endo_nbr, options_.periods);
if (options_.solve_algo == 10 || options_.solve_algo == 11)% mixed complementarity problem
    residuals(eq_to_ignore,bsxfun(@le, endogenousvariables(eq_to_ignore,:), lb(eq_to_ignore)+eps) | bsxfun(@ge,endogenousvariables(eq_to_ignore,:),ub(eq_to_ignore)-eps))=0;
end
maxerror = max(max(abs(residuals)));
success = ~check;

if ~success && options_.debug
    dprintf('solve_stacked_problem: Nonlinear solver routine failed with errorcode=%i.', errorcode)
end
