// --+ options: json=compute, stochastic +--

var y x z x1 z1 ;

varexo ey ex ez;

parameters alpha beta;

alpha = 1.0/3.0;
beta = 0.77;

model;

  [name='eq4x1']
  x1 = .5*x1(-1) + ex;

  [name='eq4z1']
  z1 = .8*z1(-1) + ez;

  [name='eq4x']
  x = x1^2/(1+x1^2);

  [name='eq4z']
  z = z1^2/(1+z1^2);

  [name='eq4y']
  y = alpha*x + (1-alpha)*z + beta*((1-alpha)*y(-1)+alpha*y(-2))+ey;

end;

/*
** Artificial dataset
*/

shocks;
  var ex = .1;
  var ez = .1;
  var ey = .1;
end;

histval;
  x1(0) = 0;
  z1(0) = 0;
  x(0) = 0;
  z(0) = 0;
  y(0) = 0;
  y(-1) = 0;
end;

simulations = simul_backward_model([], 102);

/*
** Estimation by NLS
*/

eparams.alpha = .7;
eparams.beta = 1;

simulations.ey = dseries(NaN); // Reset residuals of the equation to NaN.

if ~isoctave
    % Under Octave, estimate.nls (provided by matlab/+estimate/nls.m) is not
    % accessible because there is a function which has the same name as the
    % +estimate package, namely matlab/cli/estimate.m. This is a known Octave
    % bug (https://savannah.gnu.org/bugs/?func=detailitem&item_id=46889).
    estimate.nls('eq4y', eparams, simulations, dates('3Y'):simulations.dates(end));
end
