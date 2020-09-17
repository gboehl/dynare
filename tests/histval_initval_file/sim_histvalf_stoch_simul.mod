// Uses autonomous system from sim_base.mod, but adds separate system where exogenous variables have several leads and lags
// Lags and leads on exogenous variables are substituted out by auxiliary variables

var c cmav k z_backward z_forward;
varexo x;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;

model;
  c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
  c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
  z_backward=0.1*1+0.9*z_backward(-1) + (x(-4) - 1);
  z_forward=0.2*1+0.8*z_forward(+1) + (x(+4) - 1);
  cmav = 0.2*(c(-2) + c(-1) + c + c(+1) + c(+2));
end;

initval;
  x = 1;
end;

steady_state_model;
  k = ((bet + delt)/(aa*alph*x))^(1/(alph - 1));
  c = aa*x*k^alph - delt*k;
  z_backward = x;
  z_forward = x;
  cmav = c;
end;

steady;

shocks;
  var x;
  stderr 0.01;
end;

s = rng;
stoch_simul(periods=20, drop=0, irf=0);

reference = oo_.endo_simul;

data1 = repmat([oo_.steady_state' 1], 4, 1);
ds = dseries(data1, '1Y', [M_.endo_names; M_.exo_names]);

histval_file(series = ds);

rng(s);
stoch_simul(periods=20, drop=0, irf=0);

if max(max(abs(reference(1:5,5:end) - oo_.endo_simul(1:5,5:end)))) > 1e-8
    error('Simulation with leads and lags doesn''t match the reference')
end

data1 = repmat([oo_.steady_state' 1], 6, 1);
ds1 = dseries(data1, '1Y', [M_.endo_names; M_.exo_names]);

histval_file(series = ds1, first_obs = 6, last_obs = 6, nobs = 1);

rng(s);
stoch_simul(periods=20, drop=0, irf=0);

if max(max(abs(reference(1:5,:) - oo_.endo_simul(1:5,:)))) > 1e-8
    error('Simulation with leads and lags doesn''t match the reference')
end

histval_file(series = ds1, first_simulation_period = 7);

rng(s);
stoch_simul(periods=20, drop=0, irf=0);

if max(max(abs(reference(1:5,:) - oo_.endo_simul(1:5,:)))) > 1e-8
    error('Simulation with leads and lags doesn''t match the reference')
end

histval_file(series = ds1, first_simulation_period = 7Y);

rng(s);
stoch_simul(periods=20, drop=0, irf=0);

if max(max(abs(reference(1:5,:) - oo_.endo_simul(1:5,:)))) > 1e-8
    error('Simulation with leads and lags doesn''t match the reference')
end
