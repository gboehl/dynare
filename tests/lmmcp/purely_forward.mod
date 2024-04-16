// Regression test for bug #1720 (in the purely forward case)
// Also tests the obsolete “mcp” tag syntax

var y;

varexo eps;

model;
  [ mcp='y>1' ]
  y = sqrt(y(1))*exp(eps);
end;

initval;
  y = 1;
end;

steady;

check;


shocks;
    var eps;
    periods 1 10;
    values 1 -1;
end;

perfect_foresight_setup(periods=20);
perfect_foresight_solver(lmmcp);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

if any(oo_.endo_simul < 1)
  error('y>1 constraint not enforced')
end
