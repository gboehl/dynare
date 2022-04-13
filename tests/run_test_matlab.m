% Copyright © 2011-2017 Dynare Team
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
addpath([top_test_dir filesep 'utils']);
addpath([top_test_dir filesep '..' filesep 'matlab']);

% Test Dynare Version
if ~strcmp(dynare_version(), getenv('DYNARE_VERSION'))
  error('Incorrect version of Dynare is being tested')
end

% Test MOD files listed in Makefile.am
[modfile, name] = strtok(getenv('FILESTEM'));
[directory, testfile, ext] = fileparts([top_test_dir '/' modfile]);
cd(directory);

disp('');
disp(['***  TESTING: ' modfile ' ***']);

tic;
save(['wsMat' testfile '.mat']);
try
  dynare([testfile ext], 'console')
  testFailed = false;
catch exception
  printMakeCheckMatlabErrMsg(strtok(getenv('FILESTEM')), exception);
  testFailed = true;
end
top_test_dir = getenv('TOP_TEST_DIR');
[modfile, name] = strtok(getenv('FILESTEM'));
[directory, testfile, ext] = fileparts([top_test_dir '/' modfile]);
load(['wsMat' testfile '.mat']);
ecput = toc;
delete(['wsMat' testfile '.mat']);

cd(top_test_dir);
name = strtok(getenv('FILESTEM'));
fid = fopen([name '.m.trs'], 'w');
if fid < 0
  wd = pwd
  filestep = getenv('FILESTEM')
  error(['ERROR: problem opening file ' name '.m.trs for writing....']);
end
if testFailed
  fprintf(fid,':test-result: FAIL\n');
  fprintf(fid,':number-tests: 1\n');
  fprintf(fid,':number-failed-tests: 1\n');
  fprintf(fid,':list-of-failed-tests: %s\n', [name '.mod']);
else
  fprintf(fid,':test-result: PASS\n');
  fprintf(fid,':number-tests: 1\n');
  fprintf(fid,':number-failed-tests: 0\n');
  fprintf(fid,':list-of-passed-tests: %s\n', [name '.mod']);
end
fprintf(fid,':elapsed-time: %f\n', ecput);
fclose(fid);
warning off
exit
