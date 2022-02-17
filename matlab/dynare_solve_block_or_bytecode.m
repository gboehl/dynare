function [x,info] = dynare_solve_block_or_bytecode(y, exo, params, options, M)
% Copyright (C) 2010-2022 Dynare Team
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
        ss = x;
        if M.block_structure_stat.block(b).Simulation_Type ~= 1 && ...
                M.block_structure_stat.block(b).Simulation_Type ~= 2
            if options.solve_algo <= 4 || options.solve_algo >= 9
                [y, check] = dynare_solve('block_mfs_steadystate', ...
                                          ss(M.block_structure_stat.block(b).variable), ...
                                          options, b, ss, exo, params, T, M);
                if check ~= 0
                    %                    error(['STEADY: convergence
                    %                    problems in block ' int2str(b)])
                    info = 1;
                    return
                end
                ss(M.block_structure_stat.block(b).variable) = y;
            else
                n = length(M.block_structure_stat.block(b).variable);
                [ss, T, ~, check] = solve_one_boundary([M.fname '.static' ], ss, exo, ...
                                                       params, [], T, M.block_structure_stat.block(b).variable, n, 1, false, b, 0, options.simul.maxit, ...
                                                       options.solve_tolf, ...
                                                       options.slowc, 0, options.solve_algo, true, false, false, M, options, []);
                if check
                    info = 1;
                    return
                end
            end
        end
        % Compute endogenous if the block is of type evaluate forward/backward
        % Also update the temporary terms vector (needed for the dynare_solve case)
        [~, x, T, g1] = feval([M.fname '.static'], b, ss, exo, params, T);
    end
elseif options.bytecode
    if options.solve_algo >= 5 && options.solve_algo <= 8
        try
            x = bytecode('static', x, exo, params);
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
                [y, check] = dynare_solve('block_bytecode_mfs_steadystate', ...
                                          x(M.block_structure_stat ...
                                            .block(b).variable), ...
                                          options, b, x, exo, params, T, M);
                if check
                    %                    error(['STEADY: convergence problems in block '
                    %                    int2str(b)])
                    info = 1;
                    return
                end
                x(M.block_structure_stat.block(b).variable) = y;
            end
            % Compute endogenous if the block is of type evaluate forward/backward
            % Also update the temporary terms vector (needed for the dynare_solve case)
            try
                [~, ~, x, T] = bytecode(x, exo, params, x, 1, x, T, 'evaluate', 'static', ...
                                        ['block = ' int2str(b)]);
            catch ME
                disp(ME.message);
                info = 1;
                return
            end
        end
    else
        [x, check] = dynare_solve('bytecode_steadystate', y, ...
                                  options, exo, params);
        if check
            %            error('STEADY: convergence problems')
            info = 1;
            return
        end
    end
end
