debug = true;

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
    disp('***  TESTING: riccatiupdate.m ***');
end

if isoctave
   addpath([top_test_dir filesep '..' filesep 'mex' filesep 'octave']);
else
   addpath([top_test_dir filesep '..' filesep 'mex' filesep 'matlab']);
end

t0 = clock;

% Set the number of experiments for time measurement
N = 5000;
% Set the dimension of the problem to be solved.
r = 50;
n = 100;
tol = 1e-15;
% Set the input arguments
% P, Q: use the fact that for any real matrix A, A'*A is positive semidefinite
P = rand(n,r);
P = P'*P; 
Q = rand(n,r);
Q = Q'*Q;
K = rand(r,n);
Z = rand(n,r);
T = rand(r,r);

% 1. Update the state vairance-covariance matrix with Matlab
tElapsed1 = 0.;
tic;
for i=1:N
   Ptmp_matlab = T*(P-K*Z*P)*transpose(T)+Q;
end
tElapsed1 = toc;
disp(['Elapsed time for the Matlab Riccati update is: ' num2str(tElapsed1) ' (N=' int2str(N) ').'])

% 2. Update the state varance-covariance matrix with the mex routine
NumberOfTests = NumberOfTests+1;
tElapsed2 = 0.;
Ptmp_fortran = P;
try
   tic;
   for i=1:N
      Ptmp_fortran = riccati_update(P, T, K, Z, Q);
   end
   tElapsed2 = toc;
   disp(['Elapsed time for the Fortran Riccati update is: ' num2str(tElapsed2) ' (N=' int2str(N) ').'])
   R = norm(Ptmp_fortran-Ptmp_matlab,1);
   if (R > tol)
      testFailed = testFailed+1;
      if debug
         dprintf('The Fortran Riccati update is wrong')
      end
   end
catch
   testFailed = testFailed+1;
   if debug
      dprintf('Fortran Riccati update failed')
   end
end

% Compare the Fortran and Matlab execution time
if debug
   if tElapsed1<tElapsed2
      skipline()
      dprintf('Matlab Riccati update is %5.2f times faster than its Fortran counterpart.', tElapsed2/tElapsed1)
      skipline()
   else
      skipline()
      dprintf('Fortran Riccati update is %5.2f times faster than its Matlab counterpart.', tElapsed1/tElapsed2)
      skipline()
   end
end

% Compare results after multiple calls
N = 50;
disp(['After 1 update using the Riccati formula, the norm-1 discrepancy is ' num2str(norm(Ptmp_fortran-Ptmp_matlab,1)) '.']);
for i=2:N
   Ptmp_matlab_ini = Ptmp_matlab;
   Ptmp_fortran_ini = Ptmp_fortran;
   Ptmp_matlab = T*(Ptmp_matlab_ini-K*Z*Ptmp_matlab_ini)*transpose(T)+Q;
   Ptmp_fortran = riccati_update(Ptmp_fortran_ini, T, K, Z, Q);
end
disp(['After ' int2str(N) ' updates using the Riccati formula, the norm-1 discrepancy is ' num2str(norm(Ptmp_fortran-Ptmp_matlab,1)) '.'])

t1 = clock;

if ~debug
    cd(getenv('TOP_TEST_DIR'));
else
    dprintf('FAILED tests: %i', testFailed)
end

if  isoctave
    fid = fopen('riccatiupdate.o.trs', 'w+');
else
    fid = fopen('riccatiupdate.m.trs', 'w+');
end
if testFailed
    fprintf(fid,':test-result: FAIL\n');
    fprintf(fid,':list-of-failed-tests: riccatiupdate.m\n');
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

