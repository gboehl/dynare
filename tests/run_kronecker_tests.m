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

source_dir = getenv('source_root');
addpath([source_dir filesep 'matlab']);

dynare_config

cd kronecker

disp('**** Testing sparse_hessian_times_B_kronecker_C MEX...')
disp('')
disp('** Test 1')
info(1) = test_kron(1, num_procs);
disp('')
disp('** Test 2')
info(2) = test_kron(2, num_procs);
disp('')
disp('**** Testing A_times_B_kronecker_C MEX...')
info(3) = test_kron(3, num_procs);

num_failed_tests = sum(~info);
tests = { 'sparse1', 'sparse2', 'dense' };
failed_tests = tests(find(~info));

if num_failed_tests > 0
    fprintf('\n*** Failed tests: %s\n', failed_tests{:})
end

quit(num_failed_tests > 0)
