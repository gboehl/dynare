function [endogenousvariables, success] = sim1_purely_forward(endogenousvariables, exogenousvariables, steadystate, M_, options_)
% Performs deterministic simulation of a purely forward model

% Copyright © 2012-2023 Dynare Team
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

if ismember(options_.solve_algo, [12,14])
    [funcs, feedback_vars_idxs] = setup_time_recursive_block_simul(M_);
else
    dynamic_resid = str2func([M_.fname '.sparse.dynamic_resid']);
    dynamic_g1 = str2func([M_.fname '.sparse.dynamic_g1']);
end

function [r, J] = block_wrapper(z, feedback_vars_idx, func, y_dynamic, x, sparse_rowval, sparse_colval, sparse_colptr, T)
    % NB: do as few computations as possible inside this function, since it is
    % called a very large number of times
    y_dynamic(feedback_vars_idx) = z;
    [~, ~, r, J] = feval(func, y_dynamic, x, M_.params, steadystate, ...
                         sparse_rowval, sparse_colval, sparse_colptr, T);
end

success = true;

for it = options_.periods:-1:1
    yf = endogenousvariables(:,it+1); % Values at next period, also used as guess value for current period
    x = exogenousvariables(it,:);
    if ismember(options_.solve_algo, [12,14])
        T = NaN(M_.block_structure.dyn_tmp_nbr);
        y_dynamic = [NaN(M_.endo_nbr, 1); yf; yf];
        for blk = 1:length(M_.block_structure.block)
            sparse_rowval = M_.block_structure.block(blk).g1_sparse_rowval;
            sparse_colval = M_.block_structure.block(blk).g1_sparse_colval;
            sparse_colptr = M_.block_structure.block(blk).g1_sparse_colptr;
            if M_.block_structure.block(blk).Simulation_Type ~= 2 % Not an evaluate backward block
                [z, check, ~, ~, errorcode] = dynare_solve(@block_wrapper, y_dynamic(feedback_vars_idxs{blk}), ...
                                                           options_.simul.maxit, options_.dynatol.f, ...
                                                           options_.dynatol.x, options_, ...
                                                           feedback_vars_idxs{blk}, funcs{blk}, y_dynamic, x, sparse_rowval, sparse_colval, sparse_colptr, T);
                if check
                    success = false;
                    if options_.debug
                        dprintf('sim1_purely_forward: Nonlinear solver routine failed with errorcode=%i in block %i and period %i.', errorcode, blk, it)
                    end
                end
                y_dynamic(feedback_vars_idxs{blk}) = z;
            end
            %% Compute endogenous if the block is of type evaluate or if there are recursive variables in a solve block.
            %% Also update the temporary terms vector.
            [y_dynamic, T] = feval(funcs{blk}, y_dynamic, x, M_.params, ...
                                   steadystate, sparse_rowval, sparse_colval, ...
                                   sparse_colptr, T);
        end
        endogenousvariables(:,it) = y_dynamic(M_.endo_nbr+(1:M_.endo_nbr));
    else
        [tmp, check, ~, ~, errorcode] = dynare_solve(@dynamic_forward_model_for_simulation, yf, ...
                                                     options_.simul.maxit, options_.dynatol.f, options_.dynatol.x, ...
                                                     options_, dynamic_resid, dynamic_g1, yf, x, M_.params, steadystate, M_.dynamic_g1_sparse_rowval, M_.dynamic_g1_sparse_colval, M_.dynamic_g1_sparse_colptr);
        if check
            success = false;
            dprintf('sim1_purely_forward: Nonlinear solver routine failed with errorcode=%i in period %i.', errorcode, it)
            break
        end
        endogenousvariables(:,it) = tmp(1:M_.endo_nbr);
    end
end

end

function [r, J] = dynamic_forward_model_for_simulation(z, dynamic_resid, dynamic_g1, ylead, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr)

endo_nbr = length(z);

y = [ NaN(endo_nbr, 1); z; ylead];

[r, T_order, T] = dynamic_resid(y, x, params, steady_state);

if nargout>1
    Jacobian = dynamic_g1(y, x, params, steady_state, sparse_rowval, ...
                          sparse_colval, sparse_colptr, T_order, T);
    J = Jacobian(:, endo_nbr+(1:endo_nbr));
end

end
