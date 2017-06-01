function [ys,check] = walsh1_old_ss_steadystate(ys,exo)
global M_ 

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
    paramname = deblank(M_.param_names(ii,:));
    eval([ paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;


%% Enter model equations here

pi = thetass-1;
en = 1/3;
eR = 1/betta;
y_k = (1/alphha)*(1/betta-1+delta);
ek = en*y_k^(-1/(1-alphha));
ec = ek*(y_k-delta);
em = ec*(a/(1-a))^(-1/b)*((thetass-betta)/thetass)^(-1/b);
ey = ek*y_k;
Xss = a*ec^(1-b)*(1+(a/(1-a))^(-1/b)*((thetass-betta)/thetass)^((b-1)/b));
Psi = (1-alphha)*(ey/en)*Xss^((b-phi1)/(1-b))*a*ec^(-b)*(1-en)^eta;
n = log(en);
k = log(ek);
m = log(em);
c = log(ec);
y = log(ey);
R = log(eR);
z = 0;
u = 0;

%% end own model equations

for iter = 1:length(M_.params) %update parameters set in the file
    eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
    varname = deblank(M_.endo_names(ii,:));
    eval(['ys(' int2str(ii) ') = ' varname ';']);
end
