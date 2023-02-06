%% CODE TO SIMULATE LOGISTIC MAP
%% DECLARATION OF ENDOGENOUS VARIABLES
var x;
parameters r s;

r=3.8;
s=0;

model;
x = r*x(-1)*(1 - x(-1)) + s*x(+1);
end;

initval;
x = 0;
end;

check;

%% DETERMINISTIC SIMULATION
perfect_foresight_setup(periods = 40);
perfect_foresight_solver(stack_solve_algo=0, maxit=100);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

dsample 40;
rplot x;
