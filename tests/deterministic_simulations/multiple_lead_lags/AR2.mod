var y_backward;

parameters rho_1 rho_2;

rho_1=0.2;
rho_2=0.1;

// Equilibrium conditions
model;
y_backward=rho_1*y_backward(-1)+rho_2*y_backward(-2);
end;

// Set starting value for solver
initval;
y_backward=1;
end;
steady;

// Set initial conditions for state variables
histval;
y_backward(0)=1;
y_backward(-1)=2;
end;


// Check the Blanchard-Kahn conditions
check;

// Deterministic simulation of the model for 200 periods
perfect_foresight_setup(periods=100);
perfect_foresight_solver(tolf=1e-12);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

// Display the path of consumption and capital
rplot y_backward;

junk=zeros(1,options_.periods+M_.maximum_lag);
junk(1)=2;
junk(2)=1;

for ii=3:options_.periods+2
    junk(ii)=M_.params(strmatch('rho_1',M_.param_names,'exact'))*junk(ii-1)+M_.params(strmatch('rho_2',M_.param_names,'exact'))*junk(ii-2);
end

if max(abs(junk(M_.maximum_lag+1:end)-oo_.endo_simul(strmatch('y_backward',M_.endo_names,'exact'),:)))>1e-10
    error('Solution of purely backwards model not correct')
end
        
