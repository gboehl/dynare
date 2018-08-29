function run_all_tests()

% Copyright (C) 2018 Dynare Team
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

r = [];

r = [r; run_this_test('var-1')];
r = [r; run_this_test('var-2')];
r = [r; run_this_test('var-3')];
r = [r; run_this_test('var-4')];

r = [r; run_this_test('trend-component-1')];
r = [r; run_this_test('trend-component-2')];
r = [r; run_this_test('trend-component-3')];
r = [r; run_this_test('trend-component-4')];
r = [r; run_this_test('trend-component-5')];
r = [r; run_this_test('trend-component-6')];
r = [r; run_this_test('trend-component-7')];
r = [r; run_this_test('trend-component-9')];
r = [r; run_this_test('trend-component-10')];
r = [r; run_this_test('trend-component-11')];

print_results(r);

function o = run_this_test(folder)
try
    cd(folder)
    tstart = tic;
    clear M_ oo_ options_
    dynare example
    elapsed = toc(tstart);
    o = {folder, true, elapsed};
    system('./clean')
    cd ..
catch
    o = {folder, false, NaN};
    system('./clean')
    cd ..
end

function print_results(r)
message = sprintf('Testsuite results (PAC model):\n');
for i = 1:size(r, 1)
    if r{i,2}
        message = sprintf('%s\n%s\t\t PASS (%ss)', message, r{i,1}, num2str(r{i,3}));
    else
        message = sprintf('%s\n%s\t\t FAILED', message, r{i,1});
    end
end
disp(message)
