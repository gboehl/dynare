debug = false;

if debug
    [top_test_dir, ~, ~] = fileparts(mfilename('fullpath'));
else
    top_test_dir = getenv('TOP_TEST_DIR');
end

addpath(sprintf('%s/matlab', top_test_dir(1:end-6)))

if ~debug
    % Test Dynare Version
    if ~strcmp(dynare_version(), getenv('DYNARE_VERSION'))
        error('Incorrect version of Dynare is being tested')
    end
end

dynare_config;

NumberOfTests = 0;
testFailed = 0;

if ~debug
    skipline()
    disp('***  TESTING: nonlinearsolvers.m ***');
end

%
% TEST
%

t0 = clock;

NumberOfTests = NumberOfTests+1;

try
    dataset = dseries('simulateddata.m');

    dcontrib --model sandbox.mod --tags zpac eq:x1 --database dataset --output results --range 2023Q1:2073Q1

    if max(abs(sum(results.z.data, 2)-dataset.z(dates('2023Q1'):dates('2073Q1')).data))>1e-5

        error('Computation of dynamic contributions failed.')
    end

catch

    testFailed = testFailed+1;

end

t1 = clock;


if ~debug
    cd(getenv('TOP_TEST_DIR'));
else
    dprintf('FAILED tests: %i', testFailed)
end

if  isoctave
    fid = fopen('contribs.o.trs', 'w+');
else
    fid = fopen('contribs.m.trs', 'w+');
end
if testFailed
    fprintf(fid,':test-result: FAIL\n');
    fprintf(fid,':list-of-failed-tests: nonlinearsolvers.m\n');
else
    fprintf(fid,':test-result: PASS\n');
end
fprintf(fid,':number-tests: %i\n', NumberOfTests);
fprintf(fid,':number-failed-tests: %i\n', testFailed);
fprintf(fid,':elapsed-time: %f\n', etime(t1, t0));
fclose(fid);

if ~debug
    exit;
end
%
% END OF TEST
%
