function [x,info] = dynare_solve_block_or_bytecode(y, exo, params, options, M)

% Copyright © 2010-2023 Dynare Team
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

info = 0;
x = y;
if options.block && ~options.bytecode
    T = NaN(M.block_structure_stat.tmp_nbr, 1);
    for b = 1:length(M.block_structure_stat.block)
        fh_static = str2func(sprintf('%s.sparse.block.static_%d', M.fname, b));
        ss = x;
        if M.block_structure_stat.block(b).Simulation_Type ~= 1 && ...
                M.block_structure_stat.block(b).Simulation_Type ~= 2
            if options.solve_algo <= 4 || options.solve_algo >= 9
                [y, errorflag] = dynare_solve('block_mfs_steadystate', ...
                                              ss(M.block_structure_stat.block(b).variable), ...
                                              options.simul.maxit, options.solve_tolf, options.solve_tolx, ...
                                              options, fh_static, b, ss, exo, params, T, M);
                if errorflag
                    info = 1;
                    return
                end
                ss(M.block_structure_stat.block(b).variable) = y;
            else
                n = length(M.block_structure_stat.block(b).variable);
                [ss, T, ~, check] = solve_one_boundary(fh_static, ss, exo, ...
                                                       params, [], T, M.block_structure_stat.block(b).variable, n, 1, false, b, 0, options.simul.maxit, ...
                                                       options.solve_tolf, ...
                                                       0, options.solve_algo, true, false, false, M, options, []);
                if check
                    info = 1;
                    return
                end
            end
        end
        % Compute endogenous if the block is of type evaluate forward/backward or if there are recursive variables in a solve block.
        % Also update the temporary terms vector (needed for the dynare_solve case)
        [x, T] = fh_static(ss, exo, params, M.block_structure_stat.block(b).g1_sparse_rowval, ...
                           M.block_structure_stat.block(b).g1_sparse_colval, ...
                           M.block_structure_stat.block(b).g1_sparse_colptr, T);
    end
elseif options.bytecode
    if options.solve_algo >= 5 && options.solve_algo <= 8
        try
            if options.block
                x = bytecode('static', 'block_decomposed', x, exo, params);
            else
                x = bytecode('static', x, exo, params);
            end
        catch ME
            disp(ME.message);
            info = 1;
            return
        end
    elseif options.block
        T = NaN(M.block_structure_stat.tmp_nbr, 1);
        for b = 1:length(M.block_structure_stat.block)
            if M.block_structure_stat.block(b).Simulation_Type ~= 1 && ...
                    M.block_structure_stat.block(b).Simulation_Type ~= 2
                [y, errorflag] = dynare_solve('block_bytecode_mfs_steadystate', ...
                                              x(M.block_structure_stat.block(b).variable), ...
                                              options.simul.maxit, options.solve_tolf, options.solve_tolx, ...
                                              options, b, x, exo, params, T, M);
                if errorflag
                    info = 1;
                    return
                end
                x(M.block_structure_stat.block(b).variable) = y;
            end
            % Compute endogenous if the block is of type evaluate forward/backward or if there are recursive variables in a solve block.
            % Also update the temporary terms vector (needed for the dynare_solve case)
            try
                [~, ~, x, T] = bytecode(x, exo, params, x, 1, x, T, 'evaluate', 'static', ...
                                        'block_decomposed', ['block = ' int2str(b)]);
            catch ME
                disp(ME.message);
                info = 1;
                return
            end
        end
    else
        [x, errorflag] = dynare_solve('bytecode_steadystate', y, ...
                                      options.simul.maxit, options.solve_tolf, options.solve_tolx, ...
                                      options, exo, params);
        if errorflag
            info = 1;
            return
        end
    end
end
