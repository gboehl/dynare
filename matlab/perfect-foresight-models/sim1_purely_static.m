function [endogenousvariables, info] = sim1_purely_static(endogenousvariables, exogenousvariables, steadystate, M, options)

% Performs deterministic simulation of a purely static model

% Copyright © 2021-2023 Dynare Team
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

if ismember(options.solve_algo, [12,14])
    [funcs, feedback_vars_idxs] = setup_time_recursive_block_simul(M);
else
    dynamicmodel = str2func(sprintf('%s.%s', M.fname, 'dynamic'));
end

function [r, J] = block_wrapper(z, feedback_vars_idx, func, y_dynamic, x, sparse_rowval, sparse_colval, sparse_colptr, T)
    % NB: do as few computations as possible inside this function, since it is
    % called a very large number of times
    y_dynamic(feedback_vars_idx) = z;
    [~, ~, r, J] = feval(func, y_dynamic, x, M.params, steadystate, ...
                         sparse_rowval, sparse_colval, sparse_colptr, T);
end

info.status = true;

y = endogenousvariables(:,1);

for it = 1:options.periods
    if ismember(options.solve_algo, [12,14])
        x = exogenousvariables(it,:);
        T = NaN(M.block_structure.dyn_tmp_nbr);
        y_dynamic = [NaN(M.endo_nbr, 1); y; NaN(M.endo_nbr, 1)];
        for blk = 1:length(M.block_structure.block)
            sparse_rowval = M.block_structure.block(blk).g1_sparse_rowval;
            sparse_colval = M.block_structure.block(blk).g1_sparse_colval;
            sparse_colptr = M.block_structure.block(blk).g1_sparse_colptr;
            if M.block_structure.block(blk).Simulation_Type ~= 1  % Not an evaluate forward block
                [z, check, ~, ~, errorcode] = dynare_solve(@block_wrapper, y_dynamic(feedback_vars_idxs{blk}), ...
                                                           options.simul.maxit, options.dynatol.f, ...
                                                           options.dynatol.x, options, ...
                                                           feedback_vars_idxs{blk}, funcs{blk}, y_dynamic, x, sparse_rowval, sparse_colval, sparse_colptr, T);
                if check
                    info.status = false;
                    if options.debug
                        dprintf('sim1_purely_static: Nonlinear solver routine failed with errorcode=%i in block %i and period %i.', errorcode, blk, it)
                    end
                end
                y_dynamic(feedback_vars_idxs{blk}) = z;
            end
            %% Compute endogenous if the block is of type evaluate or if there are recursive variables in a solve block.
            %% Also update the temporary terms vector.
            [y_dynamic, T] = feval(funcs{blk}, y_dynamic, x, M.params, ...
                                   steadystate, sparse_rowval, sparse_colval, ...
                                   sparse_colptr, T);
        end
        endogenousvariables(:,it) = y_dynamic(M.endo_nbr+(1:M.endo_nbr));
    else
        [tmp, check, ~, ~, errorcode] = dynare_solve(@dynamic_static_model_for_simulation, y, ...
                                                     options.simul.maxit, options.dynatol.f, options.dynatol.x, ...
                                                     options, dynamicmodel, exogenousvariables, M.params, steadystate, it);
        if check
            info.status = false;
            if options.debug
                dprintf('sim1_purely_static: Nonlinear solver routine failed with errorcode=%i in period %i.', errorcode, it)
            end
        end
        endogenousvariables(:,it) = tmp;
    end
    y = endogenousvariables(:,it);
end

end

function [r, J] = dynamic_static_model_for_simulation(z, dynamicmodel, x, params, steady_state, it_)

% NOTE: It is assumed that all variables appear at time t in the model.

if nargout>1
    % Compute residuals and jacobian of the full dynamic model.
    [r, J] = feval(dynamicmodel, z, x, params, steady_state, it_);
    J = J(:,1:rows(J)); % Remove derivatives with respect to shocks.
else
    % Compute residuals.
    r = feval(dynamicmodel, z, x, params, steady_state, it_);
end

end
