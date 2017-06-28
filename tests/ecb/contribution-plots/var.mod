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

mydata = dseries(rand(30, 3), 1, {'ffr', 'unrate', 'cpi'});

params.p_ffr_ffr_lag_1 = 1;
params.p_ffr_ffr_lag_2 = 2;
params.p_ffr_unrate_lag_1 = 3;
params.p_ffr_cpi_lag_1 = 4;

plot_contributions('eqnum', 'ffr', mydata, params);