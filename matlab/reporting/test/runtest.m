addpath('~/builds/dynare/matlab/utilities/general')
addpath('~/builds/dynare/matlab/modules/dates/src')
addpath('~/builds/dynare/matlab/modules/dseries/src')
addpath('../src')

initialize_dates_toolbox;
initialize_dseries_toolbox;

db_a = dseries('db_a.csv');
db_q = dseries('db_q.csv');
dc_a = dseries('dc_a.csv');
dc_q = dseries('dc_q.csv');

createReport(dc_a, dc_q, db_a, db_q);