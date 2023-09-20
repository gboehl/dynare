## Copyright © 2009-2023 Dynare Team
##
## This file is part of Dynare.
##
## Dynare is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## Dynare is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

## Implementation notes:
##
## Before every call to Dynare, the contents of the workspace is saved in
## 'wsOct', and reloaded after Dynare has finished (this is necessary since
## Dynare does a 'clear -all').
## Note that we take care of clearing oo_ before saving the workspace,
## otherwise the newly created oo_ will be scratched upon loading,
## thus making the path comparison bogus.

source_dir = getenv('source_root');
addpath([source_dir filesep 'tests' filesep 'utils']);
addpath([source_dir filesep 'matlab']);

## Test block_bytecode/ls2003.mod with various combinations of
## block/bytecode/solve_algo/stack_solve_algo
failedBlock = {};
cd block_bytecode
tic;
for blockFlag = 0:1
    for storageFlag = 0:2 % 0=M-file, 1=use_dll, 2=bytecode
        default_solve_algo = 2;
        default_stack_solve_algo = 0;
        if !blockFlag && storageFlag != 2
            solve_algos = [0:4 9];
            stack_solve_algos = [0 1 6];
        elseif blockFlag && storageFlag != 2
            solve_algos = [0:4 6:9];
            stack_solve_algos = [0:4 6];
        else
            solve_algos = 0:9;
            stack_solve_algos = 0:6;
        endif

        # Workaround for strange race condition related to the static/dynamic
        # files (especially when we switch to/from use_dll)
        if exist('+ls2003_tmp')
            rmdir('+ls2003_tmp', 's');
            pause(1)
        endif

        for i = 1:length(solve_algos)
            if !blockFlag && storageFlag == 0 && (i == 1)
                ## This is the reference simulation path against which all
                ## other simulations will be tested
                try
                    old_path = path;
                    save wsOct
                    run_ls2003(blockFlag, storageFlag, solve_algos(i), default_stack_solve_algo)
                    load wsOct
                    path(old_path);
                    y_ref = oo_.endo_simul;
                    save('test.mat','y_ref');
                catch
                    load wsOct
                    path(old_path);
                    failedBlock{size(failedBlock,2)+1} = ['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(solve_algos(i)) ',' num2str(default_stack_solve_algo) ')'];
                    printMakeCheckOctaveErrMsg(['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(solve_algos(i)) ',' num2str(default_stack_solve_algo) ')'], lasterror);
                end_try_catch
            else
                try
                    old_path = path;
                    clear oo_ # Ensure that oo_.endo_simul won’t be overwritten when loading wsOct
                    save wsOct
                    run_ls2003(blockFlag, storageFlag, solve_algos(i), default_stack_solve_algo)
                    load wsOct
                    path(old_path);
                    ## Test against the reference simulation path
                    load('test.mat','y_ref');
                    diff = oo_.endo_simul - y_ref;
                    if max(max(abs(diff))) > options_.dynatol.x
                        failedBlock{size(failedBlock,2)+1} = ['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(solve_algos(i)) ',' num2str(default_stack_solve_algo) ')'];
                        differr.message = ["ERROR: simulation path differs from the reference path" ];
                        printMakeCheckOctaveErrMsg(['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(solve_algos(i)) ',' num2str(default_stack_solve_algo) ')'], differr);
                    endif
                catch
                    load wsOct
                    path(old_path);
                    failedBlock{size(failedBlock,2)+1} = ['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(solve_algos(i)) ',' num2str(default_stack_solve_algo) ')'];
                    printMakeCheckOctaveErrMsg(['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(solve_algos(i)) ',' num2str(default_stack_solve_algo) ')'], lasterror);
                end_try_catch
            endif
        endfor
        for i = 1:length(stack_solve_algos)
            try
                old_path = path;
                clear oo_ # Ensure that oo_.endo_simul won’t be overwritten when loading wsOct
                save wsOct
                run_ls2003(blockFlag, storageFlag, default_solve_algo, stack_solve_algos(i))
                load wsOct
                path(old_path);
                ## Test against the reference simulation path
                load('test.mat','y_ref');
                diff = oo_.endo_simul - y_ref;
                if max(max(abs(diff))) > options_.dynatol.x
                    failedBlock{size(failedBlock,2)+1} = ['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(default_solve_algo) ',' num2str(stack_solve_algos(i)) ')'];
                    differr.message = ["ERROR: simulation path differs from the reference path" ];
                    printMakeCheckOctaveErrMsg(['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(default_solve_algo) ',' num2str(stack_solve_algos(i)) ')'], differr);
                endif
            catch
                load wsOct
                path(old_path);
                failedBlock{size(failedBlock,2)+1} = ['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(default_solve_algo) ',' num2str(stack_solve_algos(i)) ')'];
                printMakeCheckOctaveErrMsg(['block_bytecode' filesep 'run_ls2003.m(' num2str(blockFlag) ',' num2str(storageFlag) ',' num2str(default_solve_algo) ',' num2str(stack_solve_algos(i)) ')'], lasterror);
            end_try_catch
        endfor
    endfor
endfor

delete('wsOct');

if size(failedBlock,2) > 0
    fprintf('\n*** Failed tests: %s\n', failedBlock{:})
end

fprintf('\n*** Elapsed time (in seconds): %.1f\n\n', toc);

quit(size(failedBlock,2) > 0)

## Local variables:
## mode: Octave
## End:
