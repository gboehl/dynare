// Uses autonomous system from sim_base.mod, but adds separate system where exogenous variables have several leads and lags
// Lags and leads on exogenous variables are substituted out by auxiliary variables

data1 = repmat([1.2, 1.2, 12, 1, 1, 1], 208, 1);
data1(6, 6) = 0.9; //shock to x in period 2
ds = dseries(data1, '1Y', {'c', 'cmav', 'k', 'z_backward', 'z_forward', 'x'});

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

initval_file(series = ds);
if oo_.initval_series.dates(1) ~= dates('1Y');
  error('Wrong initial date in oo_.initval_series');
end;
if oo_.initval_series{'x'}.data(6) ~= 0.9;
  error('Wrong value for x');
end;

perfect_foresight_setup(periods=200);
perfect_foresight_solver(maxit=100);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed');
end

base_results=load(['sim_exo_lead_lag' filesep 'Output' filesep 'sim_exo_lead_lag_results.mat']);
if max(max(abs(base_results.oo_.endo_simul(1:5,:) - oo_.endo_simul(1:5,:)))) > 1e-8
    error('Simulation with leads and lags doesn''t match the one with auxiliary variables')
end

data1 = repmat([1.2, 1.2, 12, 1, 1, 1], 212, 1);
data1(8, 6) = 0.9; //shock to x in period 2
ds1 = dseries(data1, '1Y', {'c', 'cmav', 'k', 'z_backward', 'z_forward', 'x'});

initval_file(series = ds1, first_obs = 3, last_obs = 210, nobs = 208);
if oo_.initval_series.dates(1) ~= dates('3Y');
  error('Wrong initial date in oo_.initval_series');
end;
if oo_.initval_series{'x'}.data(6) ~= 0.9;
  error('Wrong value for x');
end;


perfect_foresight_setup(periods=200);
perfect_foresight_solver(maxit=100);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed');
end

base_results=load(['sim_exo_lead_lag' filesep 'Output' filesep 'sim_exo_lead_lag_results.mat']);
if max(max(abs(base_results.oo_.endo_simul(1:5,:) - oo_.endo_simul(1:5,:)))) > 1e-8
    error('Simulation with leads and lags doesn''t match the one with auxiliary variables')
end

initval_file(series = ds1, first_obs = 3Y, last_obs = 210Y, nobs = 208);
if oo_.initval_series.dates(1) ~= dates('3Y');
  error('Wrong initial date in oo_.initval_series');
end;
if oo_.initval_series{'x'}.data(6) ~= 0.9;
  error('Wrong value for x');
end;


perfect_foresight_setup(periods=200);
perfect_foresight_solver(maxit=100);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed');
end

base_results=load(['sim_exo_lead_lag' filesep 'Output' filesep 'sim_exo_lead_lag_results.mat']);
if max(max(abs(base_results.oo_.endo_simul(1:5,:) - oo_.endo_simul(1:5,:)))) > 1e-8
    error('Simulation with leads and lags doesn''t match the one with auxiliary variables')
end

initval_file(series = ds1, first_simulation_period = 7);
if oo_.initval_series.dates(1) ~= dates('3Y');
  error('Wrong initial date in oo_.initval_series');
end;
if oo_.initval_series{'x'}.data(6) ~= 0.9;
  error('Wrong value for x');
end;


perfect_foresight_setup(periods=200);
perfect_foresight_solver(maxit=100);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed');
end

base_results=load(['sim_exo_lead_lag' filesep 'Output' filesep 'sim_exo_lead_lag_results.mat']);
if max(max(abs(base_results.oo_.endo_simul(1:5,:) - oo_.endo_simul(1:5,:)))) > 1e-8
    error('Simulation with leads and lags doesn''t match the one with auxiliary variables')
end

initval_file(series = ds1, first_simulation_period = 7Y);
if oo_.initval_series.dates(1) ~= dates('3Y');
  error('Wrong initial date in oo_.initval_series');
end;
if oo_.initval_series{'x'}.data(6) ~= 0.9;
  error('Wrong value for x');
end;


perfect_foresight_setup(periods=200);
perfect_foresight_solver(maxit=100);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed');
end

base_results=load(['sim_exo_lead_lag' filesep 'Output' filesep 'sim_exo_lead_lag_results.mat']);
if max(max(abs(base_results.oo_.endo_simul(1:5,:) - oo_.endo_simul(1:5,:)))) > 1e-8
    error('Simulation with leads and lags doesn''t match the one with auxiliary variables')
end

