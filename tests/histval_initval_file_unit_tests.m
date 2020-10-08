top_test_dir = getenv('TOP_TEST_DIR');
addpath([top_test_dir filesep '..' filesep 'matlab/']);
dynare_config;

cd('histval_initval_file');
num_tests = 0;
failed_tests = {};

ds = dseries(randn(10,4));

M = struct();
M.fname = '';
M.endo_nbr = 3;
M.orig_endo_nbr = 3;
M.endo_names = {'Variable_1','Variable_2','Variable_3'};
M.exo_nbr = 1;
M.exo_names = {'Variable_4'};
M.exo_det_nbr = 0;
M.orig_maximum_lag = 2;

caller = 'INITVAL';

options = struct();
options.series = ds;
ds1 = histvalf_initvalf(caller, M, options);

failed_tests = my_assert(failed_tests, all(all(ds1 == ds)), 'basic test');
num_tests = num_tests + 1;

options = struct();
options.series = ds;
options.first_obs = 2;
ds1 = histvalf_initvalf(caller, M, options);
failed_tests = my_assert(failed_tests, ds1.init == dates('2Y'), ...
                         'init test 1');
num_tests = num_tests + 1;

options = struct();
options.series = ds;
options.first_obs = 2;
options.last_obs = 9;
ds1 = histvalf_initvalf(caller, M, options);
failed_tests = my_assert(failed_tests, ds1.init == dates('2Y'), ...
                         'first_obs last_obs test 1');
failed_tests = my_assert(failed_tests, ds1.last == dates('9Y'), ...
                         'first_obs last_obs test 2');
num_tests = num_tests + 2;

options = struct();
options.series = ds;
options.last_obs = 9;
ds1 = histvalf_initvalf(caller, M, options);
failed_tests = my_assert(failed_tests, ds1.init == dates('1Y'), ...
                         'last_obs test 1');
failed_tests = my_assert(failed_tests, ds1.last == dates('9Y'), ...
                         'last_obs test 2');
num_tests = num_tests + 2;

options = struct();
options.series = ds;
options.first_obs = 2;
options.last_obs = 9;
options.nobs = 8;
ds1 = histvalf_initvalf(caller, M, options);
failed_tests = my_assert(failed_tests, ds1.init == dates('2Y'), ...
                         'first_obs, last_obs, nobs test 1');
failed_tests = my_assert(failed_tests, ds1.last == dates('9Y'), ...
                         'first_obs, last_obs, nobs test 2');
num_tests = num_tests + 2;

options = struct();
options.series = ds;
options.last_obs = 9;
options.nobs = 8;
ds1 = histvalf_initvalf(caller, M, options);
failed_tests = my_assert(failed_tests, ds1.init == dates('2Y'), ...
                         'last_obs, nobs test 1');
failed_tests = my_assert(failed_tests, ds1.last == dates('9Y'), ...
                         'last_obs, nobs test 2');
num_tests = num_tests + 2;

options = struct();
options.series = ds;
options.first_obs = 2;
options.last_obs = 9;
options.nobs = 7;

try
    ds1 = histvalf_initvalf(caller, M, options);
    error('This test didn''t catch the error')
catch me
    if ~strcmp(me.message, strcat('INITVAL_FILE: FIST_OBS, LAST_OBS and NOBS contain', ...
                                  ' inconsistent information. Use only two of these', ...
                                  ' options.'))
        failed_tests = cat(1, failed_tests, 'Wrong nobs error message' );
    end
end
num_tests = num_tests + 1;

options = struct();
options.series = ds;
options.first_obs = 0;

try
    ds1 = histvalf_initvalf(caller, M, options);
    error('This test didn''t catch the error')
catch me
    if ~strcmp(me.message, strcat(caller, '_FILE: first_obs must be a positive number'))
        failed_tests = cat(1, failed_tests, ...
                           'Wrong first period error message');
    end
end
num_tests = num_tests + 1;

options = struct();
options.series = ds;
options.last_obs = 11;

try
    ds1 = histvalf_initvalf(caller, M, options);
    error('This test didn''t catch the error')
catch me
    if ~strcmp(me.message, strcat(caller, '_FILE: last_obs = 11 is larger than the number', ...
                                  ' of observations in the dataset (10)'))
        failed_tests = cat(1, failed_tests, ...
                           'Wrong last period error message');
    end
end
num_tests = num_tests + 1;

fh = fopen('data.m', 'w');
init__ = 'INIT__ = ''1Y'';';
fprintf(fh, [init__ '\n']);
eval(init__);
names__ = 'NAMES__ = {''x'', ''y''};';
fprintf(fh, [names__ '\n']);
eval(names__);
tex__ = 'TEX__ = {''x'', ''y''};';
fprintf(fh, [tex__ '\n']);
eval(tex__);
x = randn(10, 1);
fprintf(fh, 'x = [');
fprintf(fh, '%f ', x);
fprintf(fh, '];\n');
y = randn(10, 1);
fprintf(fh, 'y = [');
fprintf(fh, '%f ', y);
fprintf(fh, '];\n');
fclose(fh);

M.endo_nbr = 1;
M.orig_endo_nbr = 1;
M.endo_names = {'y'};
M.exo_nbr = 1;
M.exo_names = {'x'};
M.exo_det_nbr = 0;

options = struct();
options.datafile = 'data.m';
series = histvalf_initvalf('INITVAL_FILE', M, options);
failed_tests = my_assert(failed_tests, series.init == dates('1Y'), ...
                         '*.m file first_obs test');
failed_tests = my_assert(failed_tests, series.nobs == 10, ...
                         '*.m file nobs test');

save('data.mat', 'INIT__', 'NAMES__', 'TEX__', 'x', 'y');
options = struct();
options.datafile = 'data.mat';
series = histvalf_initvalf('INITVAL_FILE', M, options);
failed_tests = my_assert(failed_tests, series.init == dates('1Y'), ...
                         '*.mat file first_obs test');
failed_tests = my_assert(failed_tests, series.nobs == 10, ...
                         '*.mat file nobs test');

fh = fopen('data.csv', 'w');
fprintf(fh, 'x,y\n');
for i = 1:size(x,1)
    fprintf(fh, '%f,%f\n', x(i), y(i));
end
fclose(fh);

if (ispc && ~matlab_ver_less_than('8.2')) || (~ispc && ~matlab_ver_less_than('9.0'))
    writetable(table(x,y), 'data.xlsx')
    options = struct();
    options.datafile = 'data.xlsx';
    series = histvalf_initvalf('INITVAL_FILE', M, options);
    failed_tests = my_assert(failed_tests, series.init == dates('1Y'), ...
                             '*.xlsx file first_obs test');
    failed_tests = my_assert(failed_tests, series.nobs == 10, ...
                             '*.xlsx file nobs test');
    num_tests = num_tests + 2;
end

if ispc && ~matlab_ver_less_than('8.2') % Does not work under GNU/Linux with R2020b (Unicode issue in xlsread)
    writetable(table(x,y), 'data.xls')
    options = struct();
    options.datafile = 'data.xls';
    series = histvalf_initvalf('INITVAL_FILE', M, options);
    failed_tests = my_assert(failed_tests, series.init == dates('1Y'), ...
                             '*.xls file first_obs test');
    failed_tests = my_assert(failed_tests, series.nobs == 10, ...
                             '*.xls file nobs test');
    num_tests = num_tests + 2;
end

cd(getenv('TOP_TEST_DIR'));
fid = fopen('histval_initval_file_unit_tests.m.trs', 'w+');
num_failed_tests = length(failed_tests)
if num_failed_tests > 0
  fprintf(fid,':test-result: FAIL\n');
  fprintf(fid,':number-tests: %d\n', num_tests);
  fprintf(fid,':number-failed-tests: %d\n', num_failed_tests);
  fprintf(fid,':list-of-failed-tests: %s\n', failed_tests{:});
else
  fprintf(fid,':test-result: PASS\n');
  fprintf(fid,':number-tests: %d\n', num_tests);
  fprintf(fid,':number-failed-tests: 0\n');
end
fclose(fid);
exit;
