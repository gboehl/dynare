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

% Test Dynare Version
if ~strcmp(dynare_version(), getenv('DYNARE_VERSION'))
  error('Incorrect version of Dynare is being tested')
end

% To add default directories, empty dseries objects
dynare_config

disp('');
disp(['***  TESTING: run_reporting_test_matlab.m ***']);
try
    cd([top_test_dir filesep 'reporting']);
    db_a = dseries('db_a.csv');
    db_q = dseries('db_q.csv');
    dc_a = dseries('dc_a.csv');
    dc_q = dseries('dc_q.csv');
    runDynareReport(dc_a, dc_q, db_a, db_q);
    testFailed = false;
catch
    testFailed = true;
end

cd(getenv('TOP_TEST_DIR'));
fid = fopen('run_reporting_test_matlab.m.trs', 'w+');
if testFailed
  fprintf(fid,':test-result: FAIL\n');
  fprintf(fid,':number-tests: 1\n');
  fprintf(fid,':number-failed-tests: 1\n');
  fprintf(fid,':list-of-failed-tests: run_reporting_test_matlab.m\n');
else
  fprintf(fid,':test-result: PASS\n');
  fprintf(fid,':number-tests: 1\n');
  fprintf(fid,':number-failed-tests: 0\n');
end
fprintf(fid,':elapsed-time: %f\n',0.0);
fclose(fid);
exit;
