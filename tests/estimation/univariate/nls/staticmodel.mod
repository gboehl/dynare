// --+ options: json=compute, stochastic +--

var y x z;

varexo ey ex ez;

parameters theta gamma;

gamma = .3;

model;

  [name='eqx']
  x = x(-1)^.8*exp(ex) ;

  [name='eqz']
  z = exp(ez);

  [name='eqy']
  y = x^(1-gamma)*z^gamma + ey;

end;

histval;
  x(0) = 0.1;
  z(0) = 0.5;
end;

shocks;
  var ex = .05;
  var ez = .05;
  var ey = .01;
end;

// Simulate a sample.
simulations = simul_backward_model([], 5000);

// Set residuals in the estimated equation to NaNs.
simulations.ey = dseries(NaN);

// Set initial guess for the estimated parameter.
clear('eparams')
eparams.gamma = 0.6;

// Call estimation routine.
estimate.nls('eqy', eparams, simulations, 1001Y:5000Y, 'annealing')
