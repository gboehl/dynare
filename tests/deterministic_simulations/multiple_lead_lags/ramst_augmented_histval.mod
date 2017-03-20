//Version of the ramst.mod augmented with purely backward looking AR-process that should not 
//alter any dynamics of original

// Endogenous variables: consumption and capital
var c k y_backward;

// Exogenous variable: technology level
varexo x;

// Parameters declaration and calibration
parameters alph gam delt bet aa rho_1 rho_2;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;
rho_1=0.2;
rho_2=0.1;

// Equilibrium conditions
model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1); // Resource constraint
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam); // Euler equation
y_backward=0.2*y_backward(-1)+0.1*y_backward(-2);
end;

// Set starting value for solver at steady state
initval;
x = 1;
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k;
y_backward=1;
end;
steady;

// Set initial conditions for state variables
histval;
k(0) = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
y_backward(0)=1;
y_backward(-1)=2;
end;


// Check the Blanchard-Kahn conditions
check;

// Declare a positive technological shock in period 1
shocks;
var x;
periods 1;
values 1.2;
end;

// Deterministic simulation of the model for 200 periods
simul(periods=200);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

junk=zeros(1,options_.periods+M_.maximum_lag);
junk(1)=2;
junk(2)=1;

for ii=3:options_.periods+M_.maximum_lag+M_.maximum_lead
    junk(ii)=M_.params(strmatch('rho_1',M_.param_names,'exact'))*junk(ii-1)+M_.params(strmatch('rho_2',M_.param_names,'exact'))*junk(ii-2);
end

if max(abs(junk(M_.maximum_lag+1:end)-oo_.endo_simul(strmatch('y_backward',M_.endo_names,'exact'),1:end-M_.maximum_lead)))>1e-10
    error('Solution of purely backwards model not correct')
end
        
ramst_results=load('../../ramst_results.mat');
if max(abs(ramst_results.oo_.endo_simul(strmatch('k',ramst_results.M_.endo_names,'exact'),1:end-M_.maximum_lead)-oo_.endo_simul(strmatch('k',M_.endo_names,'exact'),1:end-M_.maximum_lead)))>1e-10
    error('Solution of forward part of the model not correct')
end

if max(abs(ramst_results.oo_.endo_simul(strmatch('c',ramst_results.M_.endo_names,'exact'),2:end-M_.maximum_lead)-oo_.endo_simul(strmatch('c',M_.endo_names,'exact'),2:end-M_.maximum_lead)))>1e-10
    error('Solution of forward part of the model not correct')
end
        
