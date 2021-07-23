debug = false;

if debug
    [top_test_dir, ~, ~] = fileparts(mfilename('fullpath'));
else
    top_test_dir = getenv('TOP_TEST_DIR');
end

addpath(sprintf('%s/matlab', top_test_dir(1:end-6)))
addpath(sprintf('%s/tests/solver-test-functions', top_test_dir(1:end-6)))

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

tolf = 1e-6;
tolx = 1e-6;
maxit = 50;
factor = 10;

auxstruct = struct();

t0 = clock;

try
    NumberOfTests = NumberOfTests+1;
    x = rosenbrock();
    [x, errorflag] = block_trust_region(@rosenbrock, x, tolf, tolx, maxit, factor, false, auxstruct);
    if errorflag
        testFailed = testFailed+1;
    end
catch
    testFailed = testFailed+1;
end

try
    NumberOfTests = NumberOfTests+1;
    x = powell1();
    [x, errorflag] = block_trust_region(@powell1, x, tolf, tolx, maxit, factor, false, auxstruct);
    if errorflag
        testFailed = testFailed+1;
    end
catch
    testFailed = testFailed+1;
end

try
    NumberOfTests = NumberOfTests+1;
    x = powell2();
    [x, errorflag] = block_trust_region(@powell2, x, tolf, tolx, maxit, factor, false, auxstruct);
    if errorflag
        testFailed = testFailed+1;
    end
catch
    testFailed = testFailed+1;
end

try
    NumberOfTests = NumberOfTests+1;
    x = wood();
    [x, errorflag] = block_trust_region(@wood, x, tolf, tolx, maxit, factor, false, auxstruct);
    if errorflag
        testFailed = testFailed+1;
    end
catch
    testFailed = testFailed+1;
end

try
    NumberOfTests = NumberOfTests+1;
    % FIXME block_trust_region is diverging if x(1)<0. Note that trust_region is not finding the
    % solution for the same initial conditions. 
    x = helicalvalley();
    x(1) = 5;
    [x, errorflag] = block_trust_region(@helicalvalley, x, tolf, tolx, maxit, factor, false, auxstruct);
    if errorflag
        testFailed = testFailed+1;
    end
catch
    testFailed = testFailed+1;
end

try
    NumberOfTests = NumberOfTests+1;
    n = 10;
    x = watson(nan(n,1));
    [x, errorflag] = block_trust_region(@watson, x, tolf, tolx, maxit, factor, false, auxstruct);
    if errorflag
        testFailed = testFailed+1;
    end
catch
    testFailed = testFailed+1;
end

try
    NumberOfTests = NumberOfTests+1;
    % FIXME block_trust_region does not work for all n. 
    n = 9;
    x = chebyquad(nan(n,1));
    [x, errorflag] = block_trust_region(@chebyquad, x, tolf, tolx, maxit, factor, false, auxstruct);
    if errorflag
        testFailed = testFailed+1;
    end
catch
    testFailed = testFailed+1;
end

try
    NumberOfTests = NumberOfTests+1;
    n = 10;
    x = brown(nan(n,1));
    [x, errorflag] = block_trust_region(@brown, x, tolf, tolx, maxit, factor, false, auxstruct);
    if errorflag
        testFailed = testFailed+1;
    end
catch
    testFailed = testFailed+1;
end

try
    NumberOfTests = NumberOfTests+1;
    n = 10;
    x = discreteboundaryvalue(nan(n,1));
    [x, errorflag] = block_trust_region(@discreteboundaryvalue, x, tolf, tolx, maxit, factor, false, auxstruct);
    if errorflag
        testFailed = testFailed+1;
    end
catch
    testFailed = testFailed+1;
end

try
    NumberOfTests = NumberOfTests+1;
    n = 10;
    x = discreteintegralequation(nan(n,1));
    [x, errorflag] = block_trust_region(@discreteintegralequation, x, tolf, tolx, maxit, factor, false, auxstruct);
    if errorflag
        testFailed = testFailed+1;
    end
catch
    testFailed = testFailed+1;
end

try
    NumberOfTests = NumberOfTests+1;
    n = 10;
    x = trigonometric(nan(n,1));
    [x, errorflag] = block_trust_region(@trigonometric, x, tolf, tolx, maxit, factor, false, auxstruct);
    if errorflag
        testFailed = testFailed+1;
    end
catch
    testFailed = testFailed+1;
end

try
    NumberOfTests = NumberOfTests+1;
    n = 10;
    x = variablydimensioned(nan(n,1));
    [x, errorflag] = block_trust_region(@variablydimensioned, x, tolf, tolx, maxit, factor, false, auxstruct);
    if errorflag
        testFailed = testFailed+1;
    end
catch
    testFailed = testFailed+1;
end

try
    NumberOfTests = NumberOfTests+1;
    n = 10;
    x = broydentridiagonal(nan(n,1));
    [x, errorflag] = block_trust_region(@broydentridiagonal, x, tolf, tolx, maxit, factor, false, auxstruct);
    if errorflag
        testFailed = testFailed+1;
    end
catch
    testFailed = testFailed+1;
end

try
    NumberOfTests = NumberOfTests+1;
    n = 10;
    x = broydenbanded(nan(n,1));
    [x, errorflag] = block_trust_region(@broydenbanded, x, tolf, tolx, maxit, factor, false, auxstruct);
    if errorflag
        testFailed = testFailed+1;
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
    fid = fopen('nonlinearsolvers.o.trs', 'w+');
else
    fid = fopen('nonlinearsolvers.m.trs', 'w+');
end
if testFailed
  fprintf(fid,':test-result: FAIL\n');
else
  fprintf(fid,':test-result: PASS\n');
end
fprintf(fid,':number-tests: %i\n', NumberOfTests);
fprintf(fid,':number-failed-tests: %i\n', testFailed);
fprintf(fid,':list-of-passed-tests: nonlinearsolvers.m\n');
fprintf(fid,':elapsed-time: %f\n', etime(t1, t0));
fclose(fid);

if ~debug
    exit;
end

