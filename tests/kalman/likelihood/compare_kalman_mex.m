function [flag] = compare_kalman_mex(experience)


rng(experience.Seed);

N = experience.NumberOfSimulations;

pp = experience.Number0fObservedVariables*experience.Scale;
mm = experience.SizeOfTheStateVector*experience.Scale;
rr = experience.NumberOfStructuralShocks*experience.Scale;
measurement_error_flag = experience.MeasurementErrors;
gend = experience.NumberOfPeriods;

%% SET VARIOUS PARAMETERS
kalman_tol = 1e-12;
riccati_tol =1e-9;

%% SET THE STATE SPACE MODEL:

% I randomly choose the mm eigenvalues of the transition matrix.
TransitionEigenvalues  = rand(mm, 1)*2-1;

% I randomly choose the mm eigenvectors of the transition matrix
TransitionEigenvectors = rand(mm, mm);

% I build the transition matrix
T = TransitionEigenvectors*diag(TransitionEigenvalues)/TransitionEigenvectors;

% I randomly choose matrix R
R = randn(mm, rr);

% I randomly choose the covariance matrix of the structural innovations
E = randn(rr, 20*rr);
Q = E*transpose(E)/(20*rr);

% If needed I randomly choose the covariance matrix of the measurement errors
if measurement_error_flag == 0
    H = zeros(pp,1);
elseif measurement_error_flag == 1
    H = rand(pp,1);
elseif measurement_error_flag == 2
    E = randn(pp,20*pp);
    H = E*transpose(E)/(20*pp);
    H = (.1*eye(pp))*H*(.1*eye(pp));
else
    disp('compare_kalman_mex: unknown option!')
end

% Set the selection vector (mf)
MF = transpose(randperm(mm));
mf = MF(1:pp);

% Compute ergodic variance of the state equation
P = lyapunov_symm(T,R*Q*R',riccati_tol,1.000001, riccati_tol);

if measurement_error_flag==0
    HH = [];
elseif measurement_error_flag==1
    HH = diag(H);
elseif measurement_error_flag==2
    HH = H;
end

rng(experience.Seed*1938);

% Build datasets
Y = simul_state_space_model(T,R,Q,mf,gend*N, HH);

if isempty(HH)
    HH = 0;
end

%
% Evaluate likelihoods
%

LIK_matlab_0 = NaN(N,1);
LIK_mex_0 = NaN(N,1);

flag = 0;
Zflag = 0;
tic;
for i=1:N
    [LIK_matlab_0(i), ~] = kalman_filter(Y(:,(i-1)*gend+1:i*gend),1,gend,zeros(mm,1),P,kalman_tol,riccati_tol,0,0,T,Q,R,HH,mf,mm,pp,rr,Zflag,0,0);
end
T_matlab_0 = toc;

tic;
for i=1:N
    [LIK_mex_0(i), ~] = kalman_filter_mex(Y(:,(i-1)*gend+1:i*gend),zeros(mm,1),P,kalman_tol,riccati_tol,T,Q,R,mf,Zflag,HH);
end
T_mex_0 = toc;

dprintf('Zflag = 0: Fortran Kalman filter is %5.2f times faster than its Matlab counterpart.', T_matlab_0/T_mex_0)

if max(abs((LIK_mex_0-LIK_matlab_0)./LIK_matlab_0))>1e-6
    dprintf("Zflag = 0: discrepancy between Matlab and Fortran Kalman filter results!")
    flag = 1;
end

LIK_matlab_1 = NaN(N,1);
LIK_mex_1 = NaN(N,1);

Zflag = 1;
Z = eye(mm);
Z = Z(mf,:);
tic;
for i=1:N
    [LIK_matlab_1(i), ~] = kalman_filter(Y(:,(i-1)*gend+1:i*gend),1,gend,zeros(mm,1),P,kalman_tol,riccati_tol,0,0,T,Q,R,HH,Z,mm,pp,rr,Zflag,0,0);
end
T_matlab_1 = toc;

tic;
for i=1:N
    [LIK_mex_1(i), ~] = kalman_filter_mex(Y(:,(i-1)*gend+1:i*gend),zeros(mm,1),P,kalman_tol,riccati_tol,T,Q,R,Z,Zflag,HH);
end
T_mex_1 = toc;

dprintf('Zflag = 1: Fortran Kalman filter is %5.2f times faster than its Matlab counterpart.', T_matlab_1/T_mex_1)

if max(abs((LIK_mex_1-LIK_matlab_1)./LIK_matlab_1))>1e-6
    dprintf("Zflag = 1: discrepancy between Matlab and Fortran Kalman filter results!")
    flag = 1;
end

if max(abs((LIK_mex_0 - LIK_mex_1)./LIK_mex_1))>1e-6
    dprintf("Zflag = 1: discrepancy between results with and without the Zflag!")
    flag = 1;
end

dprintf('Cost of Zflag==1 is %5.2f%% (matlab)', 100*(T_matlab_1-T_matlab_0)/T_matlab_0)
dprintf('Cost of Zflag==1 is %5.2f%% (mex)', 100*(T_mex_1-T_mex_0)/T_mex_0)
