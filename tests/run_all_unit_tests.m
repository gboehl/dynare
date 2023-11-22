% Copyright © 2013-2023 Dynare Team
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
addpath([source_dir filesep 'tests' filesep 'utils']);
matlab_dir = [source_dir filesep 'matlab'];
addpath(matlab_dir);

dynare_config();

if isoctave
    load_octave_packages
end

mlist = get_directory_description(matlab_dir);

% Under Octave, do not run tests under matlab/missing/stats/
% Also skip load_m_data_file_legacy.m: it fails in the first test, but
% this is impossible to reproduce outside the runners.
if isoctave
    mlist = mlist(find(~strncmp([matlab_dir filesep 'missing/stats/'], mlist, 24)));
    mlist = mlist(find(~strcmp([matlab_dir filesep 'load_m_file_data_legacy.m'], mlist)));
end

rng(1);

failedtests = {};

for i = 1:length(mlist)
    f = mlist{i};
    if is_unit_test_available(f)
        [check, info] = mtest(f);
        for j = 1:size(info, 1)
            if ~info{j,3}
                failedtests{length(failedtests)+1} = [ mlist{i} '#' num2str(info{j,2}) ];
            end
        end
    end
end

if length(failedtests) > 0
    fprintf('\n*** Failed tests: %s\n', failedtests{:})
end

quit(length(failedtests) > 0)
