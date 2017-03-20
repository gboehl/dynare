var y;

varexo eps;

model;
    y = sqrt(y(1))*exp(eps);
end;

initval;
    y = 1;
    eps = 0;
end;

steady;

check;


shocks;
    var eps;
    periods 1 2;
    values 1 -1;
end;

simul(periods=5);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

expected_y = ones(1, 6);
expected_y(2) = exp(-1);
expected_y(1) = sqrt(exp(-1))*exp(1);

if max(abs(oo_.endo_simul-expected_y))>options_.dynatol.x
    error('Wrong solution!')
end
