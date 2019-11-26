var Efficiency, efficiency;

varexo EfficiencyInnovation;

parameters rho, effstar, sigma;

/*
** Calibration
*/


rho     =  0.950;
effstar =  1.000;
sigma   =  0.0001;

external_function(name=mean_preserving_spread,nargs=2);

model(use_dll);

  // Eq. n°1:
  efficiency = rho*efficiency(-1) + sigma*EfficiencyInnovation;

  // Eq. n°2:
  Efficiency = effstar*exp(efficiency-mean_preserving_spread(rho,sigma));

end;

shocks;
var EfficiencyInnovation = 1;
end;

steady_state_model;
efficiency=0;
Efficiency=effstar;
end;

steady;

options_.ep.stochastic.order = 0;
ts = extended_path([], 100, [], options_, M_, oo_);

options_.ep.stochastic.order = 1;
sts = extended_path([], 100, [], options_, M_, oo_);

// The model is backward, we do not care about future uncertainty, extended path and stochastic extended path
// should return the same results.
if max(max(abs(ts.data-sts.data)))>pi*options_.dynatol.x
   disp('Stochastic Extended Path:: Something is wrong here (potential bug in extended_path.m)!!!')
end
