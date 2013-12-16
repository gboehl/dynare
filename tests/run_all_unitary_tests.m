top_test_dir = getenv('TOP_TEST_DIR');
addpath(top_test_dir);
addpath([top_test_dir filesep '..' filesep 'matlab']);
dynare_config([], 0);

% Test Dynare Version
if ~strcmp(dynare_version(), getenv('DYNARE_VERSION'))
    error('Incorrect version of Dynare is being tested')
end

mlist = get_directory_description('../matlab');

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
fclose(fid);
exit
