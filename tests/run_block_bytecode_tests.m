% Copyright © 2011-2024 Dynare Team
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

% Implementation notes:
%
% Before every call to Dynare, the contents of the workspace is saved in
% 'wsMat.mat', and reloaded after Dynare has finished (this is necessary since
% Dynare does a 'clear -all').
% Note that we take care of clearing oo_ before saving the workspace,
% otherwise the newly created oo_ will be scratched upon loading,
% thus making the path comparison bogus.
% We also take care of clearing the 'exception' variable in all 'catch'
% block, because otherwise the next 'load wsMat' within a 'catch' block will
% overwrite the last exception.

source_dir = getenv('source_root');
addpath([source_dir filesep 'tests' filesep 'utils']);
addpath([source_dir filesep 'matlab']);

% Test block_bytecode/ls2003.mod with various combinations of
% block/bytecode/solve_algo/stack_solve_algo
failedBlock = {};
cd block_bytecode
tic;
for blockFlag = 0:1
    for storageFlag = 0:2 % 0=M-file, 1=use_dll, 2=bytecode
        default_solve_algo = 2;
        default_stack_solve_algo = 0;
        if ~blockFlag && storageFlag ~= 2
            solve_algos = [1:4 6:9];
            stack_solve_algos = [0:3 6];
        elseif blockFlag && storageFlag ~= 2
            solve_algos = [1:4 6:9];
            stack_solve_algos = [0:4 6];
        else
            solve_algos = 1:9;
            stack_solve_algos = 0:6;
        end
        if isoctave || user_has_matlab_license('optimization_toolbox')
            solve_algos = [ 0 solve_algos ];
        end

        % Workaround for strange race condition related to the static/dynamic
        % files (especially when we switch to/from use_dll)
        if isoctave && exist('+ls2003_tmp')
            rmdir('+ls2003_tmp', 's');
            pause(1)
        end

        for i = 1:length(solve_algos)
            if ~blockFlag && storageFlag == 0 && (i == 1)
                % This is the reference simulation path against which all
                % other simulations will be tested
                try
                    old_path = path;
                    save wsMat
                    run_ls2003(blockFlag, storageFlag, solve_algos(i), default_stack_solve_algo)
                    load wsMat
                    path(old_path);
                    y_ref = oo_.endo_simul;
                    save('test.mat','y_ref');
                catch exception
                    load wsMat
                    path(old_path);
                    failedBlock{size(failedBlock,2)+1} = ['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(solve_algos(i)) ',' num2str(default_stack_solve_algo) ')'];
                    printTestError(['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(solve_algos(i)) ',' num2str(default_stack_solve_algo) ')'], exception);
                    clear exception
                end
            else
                try
                    old_path = path;
                    clear oo_ % Ensure that oo_.endo_simul won’t be overwritten when loading wsMat
                    save wsMat
                    run_ls2003(blockFlag, storageFlag, solve_algos(i), default_stack_solve_algo)
                    load wsMat
                    path(old_path);
                    % Test against the reference simulation path
                    load('test.mat','y_ref');
                    diff = oo_.endo_simul - y_ref;
                    if max(max(abs(diff))) > options_.dynatol.x
                        failedBlock{size(failedBlock,2)+1} = ['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(solve_algos(i)) ',' num2str(default_stack_solve_algo) ')'];
                        if isoctave
                            exception.message = 'ERROR: simulation path differs from the reference path';
                        else
                            exception = MException('Dynare:simerr', 'ERROR: simulation path differs from the reference path');
                        end
                        printTestError(['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(solve_algos(i)) ',' num2str(default_stack_solve_algo) ')'], exception);
                        clear exception
                    end
                catch exception
                    load wsMat
                    path(old_path);
                    failedBlock{size(failedBlock,2)+1} = ['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(solve_algos(i)) ',' num2str(default_stack_solve_algo) ')'];
                    printTestError(['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(solve_algos(i)) ',' num2str(default_stack_solve_algo) ')'], exception);
                    clear exception
                end
            end
        end
        for i = 1:length(stack_solve_algos)
            try
                old_path = path;
                clear oo_ % Ensure that oo_.endo_simul won’t be overwritten when loading wsMat
                save wsMat
                run_ls2003(blockFlag, storageFlag, default_solve_algo, stack_solve_algos(i))
                load wsMat
                path(old_path);
                % Test against the reference simulation path
                load('test.mat','y_ref');
                diff = oo_.endo_simul - y_ref;
                if max(max(abs(diff))) > options_.dynatol.x
                    failedBlock{size(failedBlock,2)+1} = ['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(default_solve_algo) ',' num2str(stack_solve_algos(i)) ')'];
                    if isoctave
                        exception.message = 'ERROR: simulation path differs from the reference path';
                    else
                        exception = MException('Dynare:simerr', 'ERROR: simulation path difers from the reference path');
                    end
                    printTestError(['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(default_solve_algo) ',' num2str(stack_solve_algos(i)) ')'], exception);
                    clear exception
                end
            catch exception
                load wsMat
                path(old_path);
                failedBlock{size(failedBlock,2)+1} = ['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(default_solve_algo) ',' num2str(stack_solve_algos(i)) ')'];
                printTestError(['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(default_solve_algo) ',' num2str(stack_solve_algos(i)) ')'], exception);
                clear exception
            end
        end
    end
end

delete('wsMat.mat')

if size(failedBlock,2) > 0
    fprintf('\n*** Failed tests: %s\n', failedBlock{:})
end

fprintf('\n*** Elapsed time (in seconds): %.1f\n\n', toc);

quit(size(failedBlock,2) > 0)
