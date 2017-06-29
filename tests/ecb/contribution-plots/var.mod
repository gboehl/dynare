%
% Run using command: dynare var json=compute
%
var ffr, unrate, cpi;
varexo e_ffr, e_unrate, e_cpi;

model;

[eqnum='ffr']
     ffr = adl(ffr, 'p_ffr_ffr', 2) + adl(unrate, 'p_ffr_unrate', 1) + adl(cpi, 'p_ffr_cpi', 1);

[eqnum='unrate']
     unrate = adl(unrate, 'p_ffr_unrate', 6) + adl(cpi, 'p_unrate_cpi', 6);

[eqnum='cpi']
     cpi = adl(ffr, 'p_cpi_ffr', 6) + adl(cpi, 'p_cpi_cpi', 6);

end;

// Actual paths for the variables.
ds1 = dseries(randn(30, 3), 1, {'ffr', 'unrate', 'cpi'});

// Baseline paths for the variables.
ds0 = dseries(zeros(30, 3), 1, {'ffr', 'unrate', 'cpi'});

// Calibration of some ADL parameters, we need to fix the preprocessor in oder to allow the usual syntax.
idp = strmatch('p_ffr_ffr_lag_1', M_.param_names, 'exact');
M_.params(idp) = 1;
idp = strmatch('p_ffr_ffr_lag_2', M_.param_names, 'exact');
M_.params(idp) = 2;
idp = strmatch('p_ffr_unrate_lag_1', M_.param_names, 'exact');
M_.params(idp) = 3;
idp = strmatch('p_ffr_cpi_lag_1', M_.param_names, 'exact');
M_.params(idp) = 4;

plot_contributions('eqnum', 'ffr', ds1, ds0);