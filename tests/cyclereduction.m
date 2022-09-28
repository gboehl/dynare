debug = false;

if debug
    [top_test_dir, ~, ~] = fileparts(mfilename('fullpath'));
else
    top_test_dir = getenv('TOP_TEST_DIR');
end

addpath([top_test_dir filesep '..' filesep 'matlab']);

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
    disp('***  TESTING: cyclereduction.m ***');
end

matlab_cr_path = [top_test_dir filesep '..' filesep 'matlab' filesep 'missing' filesep 'mex' filesep 'cycle_reduction'];
addpath(matlab_cr_path);
cycle_reduction_matlab = @cycle_reduction;
rmpath(matlab_cr_path);

if isoctave
   addpath([top_test_dir filesep '..' filesep 'mex' filesep 'octave']);
else
   addpath([top_test_dir filesep '..' filesep 'mex' filesep 'matlab']);
end

t0 = clock;

% Set the dimension of the problem to be solved.
n = 500;
% Set the equation to be solved
A = eye(n);
B = diag(30*ones(n,1)); B(1,1) = 20; B(end,end) = 20; B = B - diag(10*ones(n-1,1),-1); B = B - diag(10*ones(n-1,1),1);
C = diag(15*ones(n,1)); C = C - diag(5*ones(n-1,1),-1); C = C - diag(5*ones(n-1,1),1);
cvg_tol = 1e-12;
X1 = zeros(n,n);
X2 = zeros(n,n);

% 1. Solve the equation with the Matlab cycle reduction algorithm
NumberOfTests = NumberOfTests+1;
tElapsed1 = 0.;
try
   tic; [X1,info] = cycle_reduction_matlab(C,B,A,cvg_tol,[0.]); tElapsed1 = toc;
   disp(['Elapsed time for the Matlab cycle reduction algorithm is: ' num2str(tElapsed1) ' (n=' int2str(n) ').'])
   R = C+B*X1+A*X1*X1;
   if (max(abs(R(:))) > cvg_tol)
      testFailed = testFailed+1;
      if debug
         dprintf('Matlab cycle_reduction solution is wrong')
      end
   end
catch
   testFailed = testFailed+1;
   if debug
      dprintf('Matlab cycle_reduction failed')
   end
end

% 2. Solve the equation with the Fortran cycle reduction algorithm
NumberOfTests = NumberOfTests+1;
tElapsed2 = 0.;
try
   tic; [X2,info] = cycle_reduction(C,B,A,cvg_tol,[0.]); tElapsed2 = toc;
   disp(['Elapsed time for the Fortran cycle reduction algorithm is: ' num2str(tElapsed2) ' (n=' int2str(n) ').'])
   R = C+B*X2+A*X2*X2;
   if (max(abs(R(:))) > cvg_tol)
      testFailed = testFailed+1;
      if debug
         dprintf('Fortran cycle_reduction solution is wrong')
      end
   end
catch
   testFailed = testFailed+1;
   if debug
      dprintf('Fortran cycle_reduction failed')
   end
end

% 3. Compare solutions of the Fortran and Matlab routines
NumberOfTests = NumberOfTests+1;
if (max(abs(X1(:) - X2(:))) > cvg_tol)
   testFailed = testFailed+1;
   if debug
      dprintf('Fortran and Matlab cycle reduction solutions differ');
   end
end

% Compare the Fortran and Matlab execution time
if debug
   if tElapsed1<tElapsed2
      skipline()
      dprintf('Matlab cyclic reduction is %5.2f times faster than its Fortran counterpart.', tElapsed2/tElapsed1)
      skipline()
   else
      skipline()
      dprintf('Fortran cyclic reduction is %5.2f times faster than its Matlab counterpart.', tElapsed1/tElapsed2)
      skipline()
   end
end

t1 = clock;

if ~debug
    cd(getenv('TOP_TEST_DIR'));
else
    dprintf('FAILED tests: %i', testFailed)
end

if  isoctave
    fid = fopen('cyclereduction.o.trs', 'w+');
else
    fid = fopen('cyclereduction.m.trs', 'w+');
end
if testFailed
    fprintf(fid,':test-result: FAIL\n');
else
    fprintf(fid,':test-result: PASS\n');
end
fprintf(fid,':number-tests: %i\n', NumberOfTests);
fprintf(fid,':number-failed-tests: %i\n', testFailed);
fprintf(fid,':list-of-passed-tests: cyclereduction.m\n');
fprintf(fid,':elapsed-time: %f\n', etime(t1, t0));
fclose(fid);

if ~debug
    exit;
end

