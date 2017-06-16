var y;

varexo eps;

parameters rho;

rho = 0.9;

model;
    y = y(-1)^rho*exp(eps);
end;

initval;
    y = 1;
    eps = 0;
end;

steady;

check;

shocks;
    var eps;
    periods 1;
    values 1;
end;

simul(periods=10);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

if max(abs(y-[1; exp(cumprod([1; rho*ones(9, 1)]))]))>options_.dynatol.x
    error('Wrong solution!')
end
