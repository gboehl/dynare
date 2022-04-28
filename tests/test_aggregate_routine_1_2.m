top_test_dir = getenv('TOP_TEST_DIR');
addpath([top_test_dir filesep() '..' filesep() 'matlab']);

% Test Dynare Version
if ~strcmp(dynare_version(), getenv('DYNARE_VERSION'))
  error('Incorrect version of Dynare is being tested')
end

% To add default directories, empty dseries objects
dynare_config();

disp('');
disp(['***  TESTING: test_aggregate_routine_1_2.m.trs ***']);
try
    aggregate('toto2.mod', {}, '', 'ecb/aggregate/1', 'ecb/aggregate/2');
    testFailed = false;
catch
    testFailed = true;
end

cd(getenv('TOP_TEST_DIR'));
fid = fopen('test_aggregate_routine_1_2.m.trs', 'w+');
if testFailed
  fprintf(fid,':test-result: FAIL\n');
  fprintf(fid,':number-tests: 1\n');
  fprintf(fid,':number-failed-tests: 1\n');
  fprintf(fid,':list-of-failed-tests: test_aggregate_routine_1_2.m\n');
else
  fprintf(fid,':test-result: PASS\n');
  fprintf(fid,':number-tests: 1\n');
  fprintf(fid,':number-failed-tests: 0\n');
  fprintf(fid,':list-of-passed-tests: test_aggregate_routine_1_2.m\n');
end
fprintf(fid,':elapsed-time: %f\n',0.0);
fclose(fid);
exit;