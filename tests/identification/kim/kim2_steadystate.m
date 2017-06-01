function [ys,check1]=kim2_steadystate(ys,exo)

global M_ 

for j=1:size(M_.param_names,1)
    eval([deblank(M_.param_names(j,:)),' = M_.params(j);'])
    assignin('base',deblank(M_.param_names(j,:)),M_.params(j));
end
for j=1:size(M_.endo_names,1)
    eval([deblank(M_.endo_names(j,:)),' = NaN;'])
end

check1=0;

s=betae*delta*alph/(1-betae+delta*betae);
a=as; %as^((1-alph)/(1+theta))*(delta^((phi+theta+1)/(theta+1))/s)^alph;
k=(delta/s/a)^(1/(alph-1));
i=delta*k;
c=(((a*k^alph)^(1+theta)-s*(i/s)^(1+theta))/(1-s))^(1/(1+theta))*(1-s);
lam = (1-s)^theta/c^(1+theta)/(1+theta);

%% end own model equations

for iter = 1:length(M_.params) %update parameters set in the file
    eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
    varname = deblank(M_.endo_names(ii,:));
    eval(['ys(' int2str(ii) ') = ' varname ';']);
end
