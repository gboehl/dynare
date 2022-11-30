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

% List of function handles
objfun = { @rosenbrock,
           @powell1,
           @powell2,
           @wood,
           @helicalvalley,
           @watson,
           @chebyquad,
           @brown,
           @discreteboundaryvalue,
           @discreteintegralequation,
           @trigonometric,
           @variablydimensioned,
           @broydentridiagonal,
           @broydenbanded };

% FIXME block_trust_region (mex) or trust_region (matlab) do not work for all n (not sure we can fix that).

%
% Test mex routine
%

t0 = clock;

for i=1:length(objfun)
    NumberOfTests = NumberOfTests+1;
    switch func2str(objfun{i})
         case 'helicalvalley'
           % FIXME block_trust_region is diverging if x(1)<0.
           x = helicalvalley();
           x(1) = 5;
      case 'chebyquad'
        % Fails with a system of 10 equations. 
        x = objfun{i}(nan(9,1));
      case {'watson', 'brown', 'discreteintegralequation', 'discreteboundaryvalue', 'chebyquad', 'trigonometric', 'variablydimensioned', 'broydenbanded', 'broydentridiagonal'}
        x = objfun{i}(nan(4,1));
      otherwise
        x = objfun{i}();
    end
    try
        [x, errorflag, exitflag] = block_trust_region(objfun{i}, x, tolf, tolx, maxit, factor, true, false);
        if isequal(func2str(objfun{i}), 'powell2')
            if ~errorflag
                testFailed = testFailed+1;
                if debug
                    dprintf('Nonlinear solver is expected to fail on %s function but did not return an error.', func2str(objfun{i}))
                end
            end
        elseif isequal(func2str(objfun{i}), 'trigonometric')
            % FIXME block_trust_region (mex) fails, with exit code equal to 4, but not trust_region (matlab). Would be nice to undertsand the difference. 
            if ~errorflag
                testFailed = testFailed+1;
                if debug
                    dprintf('Nonlinear solver is expected to fail on %s function but did not return an error.', func2str(objfun{i}))
                end
            end
        else
            if errorflag || norm(objfun{i}(x))>tolf
                testFailed = testFailed+1;
                if debug
                    dprintf('Nonlinear solver (mex) failed on %s function (norm(f(x))=%s).', func2str(objfun{i}), num2str(norm(objfun{i}(x))))
                end
            end
        end
    catch
        testFailed = testFailed+1;
        if debug
            dprintf('Nonlinear solver (mex) failed on %s function.', func2str(objfun{i}))
        end
    end
end

t1 = clock; etime(t1, t0)

%
% Test matlab routine
%

for i=1:length(objfun)
    NumberOfTests = NumberOfTests+1;
    switch func2str(objfun{i})
      case 'chebyquad'
        % Fails with a system of 10 equations. 
        x = objfun{i}(nan(9,1));
      case {'watson', 'brown', 'discreteintegralequation', 'discreteboundaryvalue', 'trigonometric', 'variablydimensioned', 'broydenbanded', 'broydentridiagonal'}
        x = objfun{i}(nan(10,1));
      otherwise
        x = objfun{i}();
    end
    try
        [x, errorflag, info] = trust_region(objfun{i}, x, 1:length(x), 1:length(x), true, [], tolf, tolx, maxit, factor);
        if isequal(func2str(objfun{i}), 'powell2')
            if ~errorflag
                testFailed = testFailed+1;
                if debug
                    dprintf('Nonlinear solver is expected to fail on %s function but did not return an error.', func2str(objfun{i}))
                end
            end
            if info~=3
                testFailed = testFailed+1;
                if debug
                    dprintf('Nonlinear solver is expected to fail on %s function with info==3 but did not the correct value of info.', func2str(objfun{i}))
                end
            end
        else
            if errorflag
                testFailed = testFailed+1;
                if debug
                    dprintf('Nonlinear solver failed on %s function (info=%s).', func2str(objfun{i}), int2str(info))
                end
            end
        end
    catch
        testFailed = testFailed+1;
        if debug
            dprintf('Nonlinear solver failed on %s function.', func2str(objfun{i}))
        end
    end
end

t2 = clock; etime(t2, t1)

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
    fprintf(fid,':list-of-failed-tests: nonlinearsolvers.m\n');
else
    fprintf(fid,':test-result: PASS\n');
end
fprintf(fid,':number-tests: %i\n', NumberOfTests);
fprintf(fid,':number-failed-tests: %i\n', testFailed);
fprintf(fid,':elapsed-time: %f\n', etime(t2, t0));
fclose(fid);

if ~debug
    exit;
end

