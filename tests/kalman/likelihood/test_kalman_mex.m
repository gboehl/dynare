debug = false;

if ~debug
    source_dir = getenv('source_root');
    addpath([source_dir filesep 'matlab']);
    dynare_config;
    addpath([source_dir filesep 'tests' filesep 'kalman' filesep 'likelihood']);
end

testFailed = 0;

if ~debug
    skipline()
    disp('***  TESTING: test_kalman_mex.m ***');
end

t0 = clock;

dprintf('Test 1: No measurement error')

Experience.Number0fObservedVariables = 10;
Experience.SizeOfTheStateVector = 100;
Experience.NumberOfStructuralShocks = 12;
Experience.MeasurementErrors = 0;
Experience.NumberOfPeriods = 300;

try
    flag = compare_kalman_mex(Experience);
    if (flag)
        testFailed = testFailed+1;
        if debug
            dprintf('MEX and MATLAB Kalman filters lead to different results')
        end
    end
catch
    testFailed = testFailed+1;
    if debug
        dprintf('Comparison between MEX and MATLAB Kalman filters failed')
    end
end

dprintf('Test 2: measurement error with diagonal variance-covariance matrix')
Experience.MeasurementErrors = 1;

try
    flag = compare_kalman_mex(Experience);
    if (flag)
        testFailed = testFailed+1;
        if debug
            dprintf('MEX and MATLAB Kalman filters lead to different results')
        end
    end
catch
    testFailed = testFailed+1;
    if debug
        dprintf('Comparison between MEX and MATLAB Kalman filters failed')
    end
end

dprintf('Test 3: measurement error with general variance-covariance matrix')
Experience.Number0fObservedVariables = 50;
Experience.SizeOfTheStateVector = 70;
Experience.NumberOfStructuralShocks = 50;
Experience.MeasurementErrors = 2;
Experience.NumberOfPeriods = 300;

try
    flag = compare_kalman_mex(Experience);
    if (flag)
        testFailed = testFailed+1;
        if debug
            dprintf('MEX and MATLAB Kalman filters lead to different results')
        end
    end
catch
    testFailed = testFailed+1;
    if debug
        dprintf('Comparison between MEX and MATLAB Kalman filters failed')
    end
end


if ~debug
    t1 = clock;
    fprintf('\n*** Elapsed time (in seconds): %.1f\n\n', etime(t1, t0));
    quit(testFailed > 0)
end
