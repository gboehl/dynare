function [y, success, maxerror, per_block_status] = solve_block_decomposed_problem(y, exo_simul, steady_state, options_, M_)
% Computes deterministic simulation with block option without bytecode
%
% INPUTS
%   y                [matrix]    initial path of endogenous (typically oo_.endo_simul)
%   exo_simul        [matrix]    path of exogenous
%   steady_state     [vector]    value used for the STEADY_STATE() operator
%   options_         [struct]    global options structure
%   M_               [struct]    global model structure
%
% OUTPUTS
%   y                [matrix]    computed path of endogenous
%   success          [boolean]   true in case of convergence, false otherwise
%   maxerror         [double]    ∞-norm of the residual
%   per_block_status [struct]    vector structure with per-block information about convergence

% Copyright © 2020-2024 Dynare Team
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

switch options_.stack_solve_algo
    case 0
        mthd='Sparse LU on stacked system';
    case {1,6}
        mthd='LBJ with LU solver';
    case 2
        mthd='GMRES on stacked system';
    case 3
        mthd='BiCGStab on stacked system';
    case 4
        mthd='Sparse LU solver with optimal path length on stacked system';
    case 7
        mthd='Solver from solve_algo option on stacked system';
    otherwise
        error('Unsupported stack_solve_algo value')
end
if options_.verbosity
    printline(41)
    dprintf('MODEL SIMULATION (method=%s):', mthd)
    skipline()
end

T=NaN(M_.block_structure.dyn_tmp_nbr, options_.periods+M_.maximum_lag+M_.maximum_lead);

maxerror = 0;
nblocks = length(M_.block_structure.block);
per_block_status = struct('success', cell(1, nblocks), 'error', cell(1, nblocks), 'iterations', cell(1, nblocks));

for blk = 1:nblocks
    if options_.bytecode
        fh_dynamic = @(y3n, x, params, ys, sparse_rowval, sparse_colval, sparse_colptr, T) bytecode_wrapper(y3n, x, params, ys, T, blk, M_, options_);
    else
        fh_dynamic = str2func(sprintf('%s.sparse.block.dynamic_%d', M_.fname, blk));
    end

    switch M_.block_structure.block(blk).Simulation_Type
        case {1, 2} % evaluate{Forward,Backward}
            if M_.block_structure.block(blk).Simulation_Type == 1
                range = M_.maximum_lag+1:M_.maximum_lag+options_.periods;
            else
                range = M_.maximum_lag+options_.periods:-1:M_.maximum_lag+1;
            end
            for it_ = range
                if it_ > 1 && it_ < size(y, 2)
                    y3n = reshape(y(:, it_+(-1:1)), 3*M_.endo_nbr, 1);
                elseif it_ > 1 % Purely backward model (in last period)
                    y3n = [ reshape(y(:, it_+(-1:0)), 2*M_.endo_nbr, 1); NaN(M_.endo_nbr, 1) ];
                elseif it_ < size(y, 2) % Purely forward model (in first period)
                    y3n = [ NaN(M_.endo_nbr, 1); reshape(y(:, it_+(0:1)), 2*M_.endo_nbr, 1) ];
                else % Static model
                    y3n = [ NaN(M_.endo_nbr, 1); y(:, it_); NaN(M_.endo_nbr, 1) ]
                end
                [y3n, T(:, it_)] = fh_dynamic(y3n, exo_simul(it_, :), M_.params, steady_state, ...
                                              M_.block_structure.block(blk).g1_sparse_rowval, ...
                                              M_.block_structure.block(blk).g1_sparse_colval, ...
                                              M_.block_structure.block(blk).g1_sparse_colptr, T(:, it_));
                y(:, it_) = y3n(M_.endo_nbr+(1:M_.endo_nbr));
            end
            success = true;
            maxblkerror = 0;
            iter = [];
        case {3, 4, 6, 7} % solve{Forward,Backward}{Simple,Complete}
            is_forward = M_.block_structure.block(blk).Simulation_Type == 3 || M_.block_structure.block(blk).Simulation_Type == 6;
            y_index = M_.block_structure.block(blk).variable(end-M_.block_structure.block(blk).mfs+1:end);
            [y, T, success, maxblkerror, iter] = solve_one_boundary(fh_dynamic, y, exo_simul, M_.params, steady_state, T, y_index, M_.block_structure.block(blk).NNZDerivatives, options_.periods, M_.block_structure.block(blk).is_linear, blk, M_.maximum_lag, options_.simul.maxit, options_.dynatol.f, cutoff, options_.stack_solve_algo, is_forward, true, false, M_, options_);
        case {5, 8} % solveTwoBoundaries{Simple,Complete}
            if ismember(options_.stack_solve_algo, [1 6])
                [y, T, success, maxblkerror, iter] = solve_two_boundaries_lbj(fh_dynamic, y, exo_simul, steady_state, T, blk, options_, M_);
            else
                [y, T, success, maxblkerror, iter] = solve_two_boundaries_stacked(fh_dynamic, y, exo_simul, steady_state, T, blk, cutoff, options_, M_);
            end
    end

    tmp = y(M_.block_structure.block(blk).variable, :);
    if any(isnan(tmp) | isinf(tmp))
        disp(['Inf or Nan value during the resolution of block ' num2str(blk)]);
        success = false;
    end

    per_block_status(blk).success = success;
    per_block_status(blk).error = maxblkerror;
    per_block_status(blk).iter = iter;

    maxerror = max(maxblkerror, maxerror);

    if ~success
        return
    end
end


function [y3n, T, r, g1b] = bytecode_wrapper(y3n, x, params, ys, T, blk, M_, options_)
    ypath = reshape(y3n, M_.endo_nbr, 3);
    xpath = [ NaN(1, M_.exo_nbr); x; NaN(1, M_.exo_nbr) ];
    [r, g1, ypath, T] = bytecode('evaluate', 'dynamic', 'block_decomposed', ['block=' int2str(blk) ], M_, options_, ypath, xpath, params, ys, 1, true, T);
    y3n = vec(ypath);
    if ismember(M_.block_structure.block(blk).Simulation_Type, [3, 4, 6, 7]) % solve{Forward,Backward}{Simple,Complete}
        g1b = spalloc(M_.block_structure.block(blk).mfs, M_.block_structure.block(blk).mfs, numel(g1));
    else
        g1b = spalloc(M_.block_structure.block(blk).mfs, 3*M_.block_structure.block(blk).mfs, numel(g1));
    end
    g1b(:, nonzeros(M_.block_structure.block(blk).bytecode_jacob_cols_to_sparse)) = g1(:, find(M_.block_structure.block(blk).bytecode_jacob_cols_to_sparse));
