function [y, success, maxerror, iter, per_block_status] = perfect_foresight_solver_core(M_, options_, oo_)

% Core function calling solvers for perfect foresight model
%
% INPUTS
% - M_                  [struct] contains a description of the model.
% - options_            [struct] contains various options.
% - oo_                 [struct] contains results
%
% OUTPUTS
% - y                   [double array] path for the endogenous variables (solution)
% - success             [logical] Whether a solution was found
% - maxerror            [double] contains the maximum absolute error
% - iter                [integer] Number of iterations of the underlying nonlinear solver (empty for non-iterative methods)
% - per_block_status    [struct] In the case of block decomposition, provides per-block solver status information (empty if no block decomposition)

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

if options_.lmmcp.status
    options_.stack_solve_algo=7;
    options_.solve_algo = 10;
elseif options_.stack_solve_algo==7 && options_.solve_algo == 11
    options_.lmmcp.status = 1; %Path solver
end

periods = options_.periods;

if options_.linear_approximation && ~(isequal(options_.stack_solve_algo,0) || isequal(options_.stack_solve_algo,7))
    error('perfect_foresight_solver: Option linear_approximation is only available with option stack_solve_algo equal to 0 or 7.')
end

if options_.endogenous_terminal_period && options_.stack_solve_algo ~= 0
    error('perfect_foresight_solver: option endogenous_terminal_period is only available with option stack_solve_algo equal to 0')
end

if options_.linear && (isequal(options_.stack_solve_algo, 0) || isequal(options_.stack_solve_algo, 7))
    options_.linear_approximation = true;
end

maxerror = [];
iter = [];
per_block_status = [];

if options_.block
    if M_.block_structure.time_recursive
        error('Internal error: can''t perform stacked perfect foresight simulation with time-recursive block decomposition')
    end
    if options_.bytecode
        try
            y = bytecode('dynamic', 'block_decomposed', oo_.endo_simul, oo_.exo_simul, M_.params, repmat(oo_.steady_state,1, periods+2), periods);
            success = true;
        catch ME
            disp(ME.message)
            if options_.no_homotopy
                error('Error in bytecode')
            end
            success = false;
        end
    else
        [y, success, maxerror, per_block_status] = solve_block_decomposed_problem(options_, M_, oo_);
    end
else
    if options_.bytecode
        try
            y = bytecode('dynamic', oo_.endo_simul, oo_.exo_simul, M_.params, repmat(oo_.steady_state, 1, periods+2), periods);
            success = true;
        catch ME
            disp(ME.message)
            if options_.no_homotopy
                error('Error in bytecode')
            end
            success = false;
        end
    else
        if M_.maximum_endo_lead == 0 && M_.maximum_endo_lag>0 && ~options_.lmmcp.status % Purely backward model
            [y, success] = sim1_purely_backward(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
        elseif M_.maximum_endo_lag == 0 && M_.maximum_endo_lead>0 && ~options_.lmmcp.status % Purely forward model
            [y, success] = sim1_purely_forward(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
        elseif M_.maximum_endo_lag == 0 && M_.maximum_endo_lead == 0 && ~options_.lmmcp.status % Purely static model
            [y, success] = sim1_purely_static(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
        else % General case
            switch options_.stack_solve_algo
              case 0
                if options_.linear_approximation
                    [y, success, maxerror] = sim1_linear(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, oo_.exo_steady_state, M_, options_);
                else
                    [y, success, maxerror, iter] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
                end
              case {1 6}
                if options_.linear_approximation
                    error('Invalid value of stack_solve_algo option!')
                end
                [y, success, maxerror, iter] = sim1_lbj(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
              case 7
                if options_.linear_approximation
                    if isequal(options_.solve_algo, 10) 
                        if options_.ramsey_policy && isfield(M_,'ramsey_model_constraints') && ~isempty(M_.ramsey_model_constraints)
                            warning('Due to ramsey_constraints you should not specify your model as model(linear)!')
                        elseif options_.lmmcp.status
                            warning('Due to lmmcp option, you should not specify your model as model(linear)!')
                        else
                            warning('It would be more efficient to set option solve_algo equal to 0!')
                        end
                    end
                    [y, success] = solve_stacked_linear_problem(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, oo_.exo_steady_state, M_, options_);
                else
                    [y, success, maxerror] = solve_stacked_problem(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
                end
              otherwise
                error('Invalid value of stack_solve_algo option!')
            end
        end
    end
end

% Some solvers do not compute the maximum error, so do it here if needed
if nargout > 2 && isempty(maxerror)
    if options_.bytecode
        residuals = bytecode('dynamic', 'evaluate', y, oo_.exo_simul, M_.params, oo_.steady_state, periods);
    else
        ny = size(y, 1);
        if M_.maximum_lag > 0
            y0 = y(:, M_.maximum_lag);
        else
            y0 = NaN(ny, 1);
        end
        if M_.maximum_lead > 0
            yT = y(:, M_.maximum_lag+periods+1);
        else
            yT = NaN(ny, 1);
        end
        yy = y(:,M_.maximum_lag+(1:periods));

        residuals = perfect_foresight_problem(yy(:), y0, yT, oo_.exo_simul, M_.params, oo_.steady_state, periods, M_, options_);
    end
    maxerror = max(max(abs(residuals)));
end
