source_dir = getenv('source_root');
addpath([source_dir filesep 'matlab']);

dynare_config;

testFailed = 0;

skipline()
disp('***  TESTING: cyclereduction.m ***');

matlab_cr_path = [source_dir filesep 'matlab' filesep 'missing' filesep 'mex' filesep 'cycle_reduction'];
addpath(matlab_cr_path);
cycle_reduction_matlab = @cycle_reduction;
rmpath(matlab_cr_path);

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
tElapsed1 = 0.;
try
   tic; [X1,info] = cycle_reduction_matlab(C,B,A,cvg_tol,[0.]); tElapsed1 = toc;
   disp(['Elapsed time for the Matlab cycle reduction algorithm is: ' num2str(tElapsed1) ' (n=' int2str(n) ').'])
   R = norm(C+B*X1+A*X1*X1,1);
   if (R > cvg_tol)
      testFailed = testFailed+1;
      dprintf('Matlab cycle_reduction solution is wrong')
   end
catch
   testFailed = testFailed+1;
   dprintf('Matlab cycle_reduction failed')
end

% 2. Solve the equation with the Fortran cycle reduction algorithm
tElapsed2 = 0.;
try
   tic; [X2,info] = cycle_reduction(C,B,A,cvg_tol,[0.]); tElapsed2 = toc;
   disp(['Elapsed time for the Fortran cycle reduction algorithm is: ' num2str(tElapsed2) ' (n=' int2str(n) ').'])
   R = norm(C+B*X2+A*X2*X2,1);
   if (R > cvg_tol)
      testFailed = testFailed+1;
      dprintf('Fortran cycle_reduction solution is wrong')
   end
catch
   testFailed = testFailed+1;
   dprintf('Fortran cycle_reduction failed')
end

% 3. Compare solutions of the Fortran and Matlab routines
if (norm(X1 - X2, 1) > cvg_tol)
   testFailed = testFailed+1;
   dprintf('Fortran and Matlab cycle reduction solutions differ');
end

% Compare the Fortran and Matlab execution time
if tElapsed1<tElapsed2
    skipline()
    dprintf('Matlab cyclic reduction is %5.2f times faster than its Fortran counterpart.', tElapsed2/tElapsed1)
    skipline()
else
    skipline()
    dprintf('Fortran cyclic reduction is %5.2f times faster than its Matlab counterpart.', tElapsed1/tElapsed2)
    skipline()
end

t1 = clock;

fprintf('\n*** Elapsed time (in seconds): %.1f\n\n', etime(t1, t0));

quit(testFailed > 0)
