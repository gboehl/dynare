## Copyright © 2013-2023 Dynare Team
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

source_dir = getenv('source_root');
addpath([source_dir filesep 'tests' filesep 'utils']);
addpath([source_dir filesep 'matlab']);

load_octave_packages

## To add default directories, empty dseries objects
dynare_config();

printf("\n***  TESTING:  run_reporting_test_octave.m ***\n");
try
    cd reporting
    db_a = dseries('db_a.csv');
    db_q = dseries('db_q.csv');
    dc_a = dseries('dc_a.csv');
    dc_q = dseries('dc_q.csv');
    runDynareReport(dc_a, dc_q, db_a, db_q);
    testFailed = false;
catch
    testFailed = true;
end

quit(testFailed)

## Local variables:
## mode: Octave
## End:
