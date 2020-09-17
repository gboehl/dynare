% By Willi Mutschler, September 26, 2016. Email: willi@mutschler.eu
function [ys,params,check] = RBCmodel_steadystate(ys,exo,M_,options_)
%% Step 0: initialize indicator and set options for numerical solver
check = 0;
options = optimset('Display','off','TolX',1e-12,'TolFun',1e-12);
params = M_.params;

%% Step 1: read out parameters to access them with their name
for ii = 1:M_.param_nbr
  eval([ M_.param_names{ii} ' = M_.params(' int2str(ii) ');']);
end

%% Step 2: Check parameter restrictions
if ETAc*ETAl<1 % parameter violates restriction (here it is artifical)
    check=1; %set failure indicator
    return;  %return without updating steady states
end

%% Step 3: Enter model equations here
A = 1;
RK = 1/BETTA - (1-DELTA);
K_O_N = (RK/(A*(1-ALFA)))^(-1/ALFA);
if K_O_N <= 0
    check = 1; % set failure indicator
    return;    % return without updating steady states
end
W = A*ALFA*(K_O_N)^(1-ALFA);
IV_O_N = DELTA*K_O_N;
Y_O_N = A*K_O_N^(1-ALFA);
C_O_N = Y_O_N - IV_O_N;
if C_O_N <= 0
    check = 1; % set failure indicator
    return;    % return without updating steady states
end

% The labor level
if ETAc == 1 && ETAl == 1
    N = (1-BETTA*B)*(C_O_N*(1-B))^-1*W/THETA/(1+(1-BETTA*B)*(C_O_N*(1-B))^-1*W/THETA);
else
    % No closed-form solution use a fixed-point algorithm
    N0 = 1/3;
    [N,~,exitflag] = fsolve(@(N) THETA*(1-N)^(-ETAl)*N^ETAc - (1-BETTA*B)*(C_O_N*(1-B))^(-ETAc)*W, N0,options);
    if exitflag <= 0
        check = 1; % set failure indicator
        return     % return without updating steady states
    end
end

C=C_O_N*N;
Y=Y_O_N*N;
IV=IV_O_N*N;
K=K_O_N*N;
LA = (C-B*C)^(-ETAc)-BETTA*B*(C-B*C)^(-ETAc);

k=log(K);
c=log(C);
a=log(A);
iv=log(IV);
y=log(Y);
la=log(LA);
n=log(N);
rk=log(RK);
w=log(W);
%% Step 4: Update parameters and variables
params=NaN(M_.param_nbr,1);
for iter = 1:M_.param_nbr %update parameters set in the file
  eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

for ii = 1:M_.orig_endo_nbr %auxiliary variables are set automatically
  eval(['ys(' int2str(ii) ') = ' M_.endo_names{ii} ';']);
end

end
