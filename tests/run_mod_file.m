% Copyright © 2011-2023 Dynare Team
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

source_root = getenv('source_root');
addpath([source_root filesep 'tests' filesep 'utils']);
addpath([source_root filesep 'matlab']);

if isoctave
    load_octave_packages
end

fprintf('\n*** TESTING: %s ***\n\n', getenv('mod_file'));

tic;

% NB: all variables will be cleared by the call to Dynare
try
    dynare(getenv('mod_file'), 'console')
    testFailed = false;
catch exception
    printTestError(getenv('mod_file'), exception);
    testFailed = true;
end


fprintf('\n*** Elapsed time (in seconds): %.1f\n\n', toc);

% Ensure proper termination and exit code
quit(testFailed)
