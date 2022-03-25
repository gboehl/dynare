function [oo_, maxerror] = perfect_foresight_solver_core(M_, options_, oo_)

% Core function calling solvers for perfect foresight model
%
% INPUTS
% - M_                  [struct] contains a description of the model.
% - options_            [struct] contains various options.
% - oo_                 [struct] contains results
%
% OUTPUTS
% - oo_                 [struct] contains results
% - maxerror            [double] contains the maximum absolute error

% Copyright (C) 2015-2022 Dynare Team
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

if options_.slowc ~= 1 && (options_.block || options_.bytecode)
    % The code is buggy and leads to wrong results, so forbid this combination
    error('Changing the value of the slowc option is not supported with block and/or bytecode option')
end

if options_.block
    if options_.bytecode
        try
            oo_.endo_simul = bytecode('dynamic', oo_.endo_simul, oo_.exo_simul, M_.params, repmat(oo_.steady_state,1, periods+2), periods);
            oo_.deterministic_simulation.status = true;
        catch ME
            disp(ME.message)
            if options_.no_homotopy
                error('Error in bytecode')
            end
            oo_.deterministic_simulation.status = false;
        end
    else
        oo_ = solve_block_decomposed_problem(options_, M_, oo_);
    end
else
    if options_.bytecode
        try
            oo_.endo_simul = bytecode('dynamic', oo_.endo_simul, oo_.exo_simul, M_.params, repmat(oo_.steady_state, 1, periods+2), periods);
            oo_.deterministic_simulation.status = true;
        catch ME
            disp(ME.message)
            if options_.no_homotopy
                error('Error in bytecode')
            end
            oo_.deterministic_simulation.status = false;
        end
    else
        if M_.maximum_endo_lead == 0 && M_.maximum_endo_lag>0 && ~options_.lmmcp.status % Purely backward model
            [oo_.endo_simul, oo_.deterministic_simulation] = ...
                sim1_purely_backward(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
        elseif M_.maximum_endo_lag == 0 && M_.maximum_endo_lead>0 && ~options_.lmmcp.status % Purely forward model
            [oo_.endo_simul, oo_.deterministic_simulation] = ...
                sim1_purely_forward(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
        elseif M_.maximum_endo_lag == 0 && M_.maximum_endo_lead == 0 && ~options_.lmmcp.status % Purely static model
            [oo_.endo_simul, oo_.deterministic_simulation] = ...
                sim1_purely_static(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
        else % General case
            switch options_.stack_solve_algo
              case 0
                if options_.linear_approximation
                    [oo_.endo_simul, oo_.deterministic_simulation] = ...
                        sim1_linear(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, oo_.exo_steady_state, M_, options_);
                else
                    [oo_.endo_simul, oo_.deterministic_simulation] = ...
                        sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
                end
              case 6
                if options_.linear_approximation
                    error('Invalid value of stack_solve_algo option!')
                end
                [oo_.endo_simul, oo_.deterministic_simulation] = ...
                    sim1_lbj(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
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
                    [oo_.endo_simul, oo_.deterministic_simulation] = ...
                        solve_stacked_linear_problem(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, oo_.exo_steady_state, M_, options_);
                else
                    [oo_.endo_simul, oo_.deterministic_simulation] = ...
                        solve_stacked_problem(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);
                end
              otherwise
                error('Invalid value of stack_solve_algo option!')
            end
        end
    end
end

if nargout>1
    if options_.lmmcp.status
        maxerror = NaN; % Could be improved
    elseif options_.block && ~options_.bytecode
        maxerror = oo_.deterministic_simulation.error;
    else
        if options_.bytecode
            residuals = bytecode('dynamic','evaluate', oo_.endo_simul, oo_.exo_simul, M_.params, oo_.steady_state, 1);
        else
            if M_.maximum_lag > 0
                y0 = oo_.endo_simul(:, M_.maximum_lag);
            else
                y0 = NaN(ny, 1);
            end
            if M_.maximum_lead > 0
                yT = oo_.endo_simul(:, M_.maximum_lag+periods+1);
            else
                yT = NaN(ny, 1);
            end
            yy = oo_.endo_simul(:,M_.maximum_lag+(1:periods));

            residuals = perfect_foresight_problem(yy(:), y0, yT, oo_.exo_simul, M_.params, oo_.steady_state, periods, M_, options_);
        end
        maxerror = max(max(abs(residuals)));
    end
end
