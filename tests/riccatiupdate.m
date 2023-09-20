debug = true;

source_dir = getenv('source_root');
addpath([source_dir filesep 'matlab']);

dynare_config;

testFailed = 0;

if ~debug
    skipline()
    disp('***  TESTING: riccatiupdate.m ***');
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
% Computing an upperbound for the norm the updated variance-covariance matrix
ub = norm(T,1)^2*norm(P,1)*(1+norm(K*Z,1))+norm(Q,1);
% Weighting the P and Q matrices to keep the norm of the variance-covariance matrix below 1
P = 0.5*P/ub;
Q = 0.5*Q/ub;

% 1. Update the state vairance-covariance matrix with Matlab
tElapsed1 = 0.;
tic;
for i=1:N
   Ptmp_matlab = T*(P-K*Z*P)*transpose(T)+Q;
end
tElapsed1 = toc;
disp(['Elapsed time for the Matlab Riccati update is: ' num2str(tElapsed1) ' (N=' int2str(N) ').'])

% 2. Update the state varance-covariance matrix with the mex routine
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
   Ptmp_matlab = T*(Ptmp_matlab-K*Z*Ptmp_matlab)*transpose(T)+Q;
   Ptmp_fortran = riccati_update(Ptmp_fortran, T, K, Z, Q);
   disp(['After ' int2str(i) ' updates using the Riccati formula, the norm-1 discrepancy is ' num2str(norm(Ptmp_fortran-Ptmp_matlab,1)) '.'])
end

t1 = clock;

fprintf('\n*** Elapsed time (in seconds): %.1f\n\n', etime(t1, t0));

quit(testFailed > 0)
