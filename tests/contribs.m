debug = false;

if ~debug
    source_dir = getenv('source_root');
    addpath([source_dir filesep 'matlab']);
end

dynare_config;

testFailed = 0;

if ~debug
    skipline()
    disp('***  TESTING: contribs.m ***');
end

%
% TEST
%

t0 = clock;

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
    fprintf('\n*** Elapsed time (in seconds): %.1f\n\n', etime(t1, t0));
    quit(testFailed > 0)
end

%
% END OF TEST
%
