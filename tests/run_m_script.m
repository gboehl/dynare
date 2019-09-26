% Copyright (C) 2015 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

top_test_dir = getenv('TOP_TEST_DIR');
[mfile, name] = strtok(getenv('FILESTEM'));

[directory, mscript, ext] = fileparts([top_test_dir '/' mfile]);
cd(directory);

try
  mscript;
  testFailed = false;
catch exception
  printMakeCheckMatlabErrMsg(strtok(getenv('FILESTEM')), exception);
  testFailed = true;
end

cd(top_test_dir);
name = strtok(getenv('FILESTEM'));
fid = fopen([name '.m.tls'], 'w');
if fid < 0
  wd = pwd
  filestep = getenv('FILESTEM')
  error(['ERROR: problem opening file ' name '.m.tls for writing....']);
end
if testFailed
  fprintf(fid,':test-result: FAIL\n');
  fprintf(fid,':number-tests: 1\n');
  fprintf(fid,':number-failed-tests: 1\n');
  fprintf(fid,':list-of-failed-tests: %s\n', [name '.m']);
else
  fprintf(fid,':test-result: PASS\n');
  fprintf(fid,':number-tests: 1\n');
  fprintf(fid,':number-failed-tests: 0\n');
  fprintf(fid,':list-of-passed-tests: %s\n', [name '.m']);
end
fclose(fid);
exit;