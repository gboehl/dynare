var y i pi ;

varexo e;

parameters a1 a2 a3 a4 a5;

a1 = -.5;
a2 =  .1;
a3 =  .9;
a4 = 1.5;
a5 = 0.5;

model;

y = y(1)+a1*(i-pi(1))+e;

pi = a2*y + a3*pi(1);

i = max(0, a4*pi+a5*y);

end;

steady_state_model;
y=0;
i=0;
pi=0;
end;

steady;
check;

shocks;
var e;
periods 1 2;
values  .3  -0.1;
end;

simul(periods=5);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

% Initialize the analytical solution for the endogenous variables.
expected_y = zeros(1, 6);
expected_pi = zeros(1, 6);
expected_i = zeros(1, 6);

% Set period 2.
tmp = inv(eye(3)-[0 0 a1; a2 0 0; a5 a4 0])*[oo_.exo_simul(2); 0; 0];

if tmp(3)>=0
    expected_y(2) = tmp(1);
    expected_pi(2) = tmp(2);
    expected_i(2) = tmp(3);
else
    expected_y(2) = oo_.exo_simul(2);
    expected_pi(2) = expected_y(2)*a2;
    expected_i(2) = 0;
end

% Set period 1.
tmp = inv(eye(3)-[0 0 a1; a2 0 0; a5 a4 0])*[oo_.exo_simul(1)+expected_y(2)-a1*expected_pi(2); a3*expected_pi(2); 0];

if tmp(3)>=0
    expected_y(1) = tmp(1);
    expected_pi(1) = tmp(2);
    expected_i(1) = tmp(3);
else
    expected_y(1) = oo_.exo_simul(1)+expected_y(2)-a1*expected_pi(2);
    expected_pi(1) = expected_y(2)*a2+a3*expected_pi(2);
    expected_i(1) = 0;
end

% Compare the paths returned by sim1_purely_forward routine and the analytical solution.
if max(abs(oo_.endo_simul(1,:)-expected_y))>options_.dynatol.x || ...
     max(abs(oo_.endo_simul(2,:)-expected_i))>options_.dynatol.x || ...
     max(abs(oo_.endo_simul(3,:)-expected_pi))>options_.dynatol.x
    error('Wrong solution!')
end
