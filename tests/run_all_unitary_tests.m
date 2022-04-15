% Copyright © 2013-2022 Dynare Team
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

top_test_dir = getenv('TOP_TEST_DIR');
addpath([top_test_dir filesep '..' filesep 'matlab']);
dynare_config();

if isoctave
    load_octave_packages
end

% Test Dynare Version
if ~strcmp(dynare_version(), getenv('DYNARE_VERSION'))
    error('Incorrect version of Dynare is being tested')
end

mlist = get_directory_description('../matlab');

% Under Octave, do not run tests under matlab/missing/stats/
% Also skip load_m_data_file_legacy.m: it fails in the first test, but
% this is impossible to reproduce outside the runners.
if isoctave
    mlist = mlist(find(~strncmp('../matlab/missing/stats/', mlist, 24)));
    mlist = mlist(find(~strcmp('../matlab/load_m_file_data_legacy.m', mlist)));
end

% Set random seed, for reproducibility
if isoctave && octave_ver_less_than('7')
    randn('state',1);
    rand('state',1);
else
    rng(1);
end

failedtests = {};

counter = 0;

for i = 1:length(mlist)
    f = [top_test_dir filesep mlist{i} ];
    if is_unitary_test_available(f)
        [check, info] = mtest(f);
        for j = 1:size(info, 1)
            counter = counter + 1;
            if ~info{j,3}
                failedtests{length(failedtests)+1} = [ mlist{i} '#' num2str(info{j,2}) ];
            end
        end
    end
end

cd(getenv('TOP_TEST_DIR'));
if isoctave
    fid = fopen('run_all_unitary_tests.o.trs', 'w+');
else
    fid = fopen('run_all_unitary_tests.m.trs', 'w+');
end
if length(failedtests) > 0
  fprintf(fid,':test-result: FAIL\n');
  fprintf(fid,':number-tests: %d\n', counter);
  fprintf(fid,':number-failed-tests: %d\n', length(failedtests));
  fprintf(fid,':list-of-failed-tests: %s\n', failedtests{:});
else
  fprintf(fid,':test-result: PASS\n');
  fprintf(fid,':number-tests: %d\n', counter);
  fprintf(fid,':number-failed-tests: 0\n');
end
fprintf(fid,':elapsed-time: %f\n',0.0);
fclose(fid);
if ~isoctave
    exit
end
