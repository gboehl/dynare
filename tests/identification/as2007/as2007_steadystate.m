function [ys,check1]=as2007_steadystate(ys,exo)

global M_

for j=1:size(M_.param_names,1)
    eval([deblank(M_.param_names(j,:)),' = M_.params(j);'])
    assignin('base',deblank(M_.param_names(j,:)),M_.params(j));
end
for j=1:size(M_.endo_names,1)
    eval([deblank(M_.endo_names(j,:)),' = NaN;'])
end

check1=0;

pie=0;
y=0;
R=0;
g=0;
z=0;
YGR=gam_steady;
INFL = pi_steady;
INT = pi_steady+rr_steady+4*gam_steady;

%% end own model equations

for iter = 1:length(M_.params) %update parameters set in the file
    eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
    varname = deblank(M_.endo_names(ii,:));
    eval(['ys(' int2str(ii) ') = ' varname ';']);
end
