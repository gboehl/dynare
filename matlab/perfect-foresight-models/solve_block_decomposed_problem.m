function oo_ = solve_block_decomposed_problem(options_, M_, oo_)
% Computes deterministic simulation with block option without bytecode

% Copyright (C) 2020-2022 Dynare Team
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

cutoff = 1e-15;

if options_.stack_solve_algo==0
    mthd='Sparse LU';
elseif options_.stack_solve_algo==1
    mthd='Relaxation';
elseif options_.stack_solve_algo==2
    mthd='GMRES';
elseif options_.stack_solve_algo==3
    mthd='BICGSTAB';
elseif options_.stack_solve_algo==4
    mthd='OPTIMPATH';
else
    mthd='UNKNOWN';
end
if options_.verbosity
    printline(41)
    disp(sprintf('MODEL SIMULATION (method=%s):',mthd))
    skipline()
end

y=oo_.endo_simul;
T=NaN(M_.block_structure.dyn_tmp_nbr, options_.periods+M_.maximum_lag+M_.maximum_lead);
oo_.deterministic_simulation.status = 0;

funcname = [ M_.fname '.dynamic'];

for blk = 1:length(M_.block_structure.block)
    recursive_size = M_.block_structure.block(blk).endo_nbr - M_.block_structure.block(blk).mfs;
    y_index = M_.block_structure.block(blk).variable((recursive_size+1):end);

    if M_.block_structure.block(blk).Simulation_Type == 1 || ... % evaluateForward
       M_.block_structure.block(blk).Simulation_Type == 2        % evaluateBackward
        oo_.deterministic_simulation.status = true;
        oo_.deterministic_simulation.error = 0;
        oo_.deterministic_simulation.iterations = 0;
        oo_.deterministic_simulation.block(blk).status = true;
        oo_.deterministic_simulation.block(blk).error = 0;
        oo_.deterministic_simulation.block(blk).iterations = 0;
        if M_.block_structure.block(blk).Simulation_Type == 1
            range = M_.maximum_lag+1:M_.maximum_lag+options_.periods;
        else
            range = M_.maximum_lag+options_.periods:-1:M_.maximum_lag+1;
        end
        for it_ = range
            y2 = dynvars_from_endo_simul(y, it_, M_);
            [~, y2, T(:, it_)] = feval(funcname, blk, y2, oo_.exo_simul, M_.params, oo_.steady_state, T(:, it_), it_, false);
            y(find(M_.lead_lag_incidence(M_.maximum_endo_lag+1, :)), it_) = y2(nonzeros(M_.lead_lag_incidence(M_.maximum_endo_lag+1, :)));
        end
    elseif M_.block_structure.block(blk).Simulation_Type == 3 || ... % solveForwardSimple
           M_.block_structure.block(blk).Simulation_Type == 4 || ... % solveBackwardSimple
           M_.block_structure.block(blk).Simulation_Type == 6 || ... % solveForwardComplete
           M_.block_structure.block(blk).Simulation_Type == 7        % solveBackwardComplete
        is_forward = M_.block_structure.block(blk).Simulation_Type == 3 || M_.block_structure.block(blk).Simulation_Type == 6;
        [y, T, oo_] = solve_one_boundary(funcname, y, oo_.exo_simul, M_.params, oo_.steady_state, T, y_index, M_.block_structure.block(blk).NNZDerivatives, options_.periods, M_.block_structure.block(blk).is_linear, blk, M_.maximum_lag, options_.simul.maxit, options_.dynatol.f, options_.slowc, cutoff, options_.stack_solve_algo, is_forward, true, false, M_, options_, oo_);
    elseif M_.block_structure.block(blk).Simulation_Type == 5 || ... % solveTwoBoundariesSimple
           M_.block_structure.block(blk).Simulation_Type == 8        % solveTwoBoundariesComplete
        [y, T, oo_] = solve_two_boundaries(funcname, y, oo_.exo_simul, M_.params, oo_.steady_state, T, y_index, M_.block_structure.block(blk).NNZDerivatives, options_.periods, M_.block_structure.block(blk).maximum_lag, M_.block_structure.block(blk).maximum_lead, M_.block_structure.block(blk).is_linear, blk, M_.maximum_lag, options_.simul.maxit, options_.dynatol.f, options_.slowc, cutoff, options_.stack_solve_algo, options_, M_, oo_);
    end

    tmp = y(M_.block_structure.block(blk).variable, :);
    if any(isnan(tmp) | isinf(tmp))
        disp(['Inf or Nan value during the resolution of block ' num2str(blk)]);
        oo_.deterministic_simulation.status = false;
        oo_.deterministic_simulation.error = 100;
        oo_.deterministic_simulation.block(blk).status = false;
        oo_.deterministic_simulation.block(blk).error = 100;
    end
    if ~oo_.deterministic_simulation.status
        return
    end
end

oo_.endo_simul = y;
