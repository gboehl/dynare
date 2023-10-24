function [flag] = compare_kalman_mex(experience)

pp = experience.Number0fObservedVariables;
mm = experience.SizeOfTheStateVector;
rr = experience.NumberOfStructuralShocks;
measurement_error_flag = experience.MeasurementErrors;
gend = experience.NumberOfPeriods;

%% SET VARIOUS PARAMETERS 
kalman_tol = 1e-12;
riccati_tol =1e-9;

%% SET THE STATE SPACE MODEL:

% I randomly choose the mm eigenvalues of the transition matrix.
TransitionEigenvalues  = rand(mm,1)*2-1;
% I randomly choose the mm eigenvectors of the transition matrix
tmp = rand(mm,mm*100);
TransitionEigenvectors = tmp*tmp'/(mm*100);
TransitionEigenvectors = rand(mm,mm);
% I build the transition matrix
T = TransitionEigenvectors*diag(TransitionEigenvalues)/TransitionEigenvectors;
% I randomly choose matrix R
R = randn(mm,rr);
% I randomly choose the covariance matrix of the structural innovations
E = randn(rr,20*rr);
Q = E*transpose(E)/(20*rr);
% If needed I randomly choose the covariance matrix of the measurement errors 
if measurement_error_flag == 0
    H = zeros(pp,1);
elseif measurement_error_flag == 1
    H = rand(pp,1);
elseif measurement_error_flag == 2
    E = randn(pp,20*pp);
    H = E*transpose(E)/(20*pp);
else
    disp('compare_kalman_mex: unknown option!')
end
% Set the selection vector (mf) 
MF = transpose(randperm(mm));
mf = MF(1:pp);

P = lyapunov_symm(T,R*Q*R',riccati_tol,1.000001, riccati_tol);

%% BUILD DATA SET (zero mean):
a = zeros(mm,1);
if measurement_error_flag == 0
    Y = simul_state_space_model(T,R,Q,mf,gend);
elseif measurement_error_flag == 1
    H = rand(pp,1);
    Y = simul_state_space_model(T,R,Q,mf,gend,diag(H));
elseif measurement_error_flag == 2
    E = randn(pp,20*pp);
    H = E*transpose(E)/(20*pp);
    Y = simul_state_space_model(T,R,Q,mf,gend,H);
else
    disp('compare_kalman_mex: unknown option!');
end

if measurement_error_flag==0
   HH = 0;
elseif measurement_error_flag==1
   HH = diag(H);
elseif measurement_error_flag==2
   HH = H;
end

flag = 0;
Zflag = 0;
tic;
[LIK_matlab,lik_matlab] = kalman_filter(Y,1,gend,a,P,kalman_tol,riccati_tol,0,0,T,Q,R,HH,mf,mm,pp,rr,Zflag,0,0);
T_matlab = toc;

tic;
[LIK_mex,lik_mex] = kalman_filter_mex(Y,a,P,kalman_tol,riccati_tol,T,Q,R,mf,Zflag,HH);
T_mex = toc;

if T_matlab<T_mex
   dprintf('Zflag = 0: Matlab Kalman filter is %5.2f times faster than its Fortran counterpart.', T_mex/T_matlab)
else
   dprintf('Zflag = 0: Fortran Kalman filter is %5.2f times faster than its Matlab counterpart.', T_matlab/T_mex)
end

if ((abs((LIK_matlab - LIK_mex)/LIK_matlab) > 1e-6) || (max(abs((lik_matlab-lik_mex)./lik_matlab)) > 1e-6))
   dprintf("Zflag = 0: discrepancy between Matlab and Fortran Kalman filter results!")
   flag = 1;
end

Zflag = 1;
Z = eye(mm);
Z = Z(mf,:);
tic;
[LIK_matlab_z,lik_matlab_z] = kalman_filter(Y,1,gend,a,P,kalman_tol,riccati_tol,0,0,T,Q,R,HH,Z,mm,pp,rr,Zflag,0,0);
T_matlab = toc;

tic;
[LIK_mex_z,lik_mex_z] = kalman_filter_mex(Y,a,P,kalman_tol,riccati_tol,T,Q,R,Z,Zflag,HH);
T_mex = toc;

if T_matlab<T_mex
   dprintf('Zflag = 1: Matlab Kalman filter is %5.2f times faster than its Fortran counterpart.', T_mex/T_matlab)
else
   dprintf('Zflag = 1: Fortran Kalman filter is %5.2f times faster than its Matlab counterpart.', T_matlab/T_mex)
end

if ((abs((LIK_matlab_z - LIK_mex_z)/LIK_matlab_z) > 1e-6) || (max(abs((lik_matlab_z-lik_mex_z)./lik_matlab_z)) > 1e-6))
   dprintf("Zflag = 1: discrepancy between Matlab and Fortran Kalman filter results!")
   flag = 1;
end

if ((abs((LIK_mex - LIK_mex_z)/LIK_mex) > 1e-6) || (max(abs((lik_mex-lik_mex_z)./lik_mex)) > 1e-6))
   dprintf("Zflag = 1: discrepancy between results with and without the Zflag!")
   flag = 1;
end
