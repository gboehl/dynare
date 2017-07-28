// --+ options: json=compute +--

/* REMARK
** ------
**
** You need to have the first line on top of the mod file. The options defined on this line are passed
** to the dynare command (you can add other options, separated by spaces or commas). The option defined
** here is mandatory for the decomposition. It forces Dynare to output another representation of the
** model in a json file (additionaly to the matlab files) which is used here to manipulate the equation
** of interest.
*/

// Declaration of the endogenous variables
var ffr, unrate, cpi;

// Declaration of the exogenous variables
varexo e_ffr, e_unrate, e_cpi;


// Declaration of the model. Note that you must associate a name to the equations. This is mandatory to
// select the equation for which you need to perform the decomposition.
model;

[name='ffr']
     ffr = adl(ffr, 'p_ffr_ffr', [1:3]) + adl(unrate, 'p_ffr_unrate', 1) + adl(cpi, 'p_ffr_cpi', [4]);

[name='unrate']
     unrate = adl(unrate, 'p_ffr_unrate', [4 2 5]) + adl(cpi, 'p_unrate_cpi', 6);

[name='cpi']
     cpi = adl(ffr, 'p_cpi_ffr', 2) + adl(cpi, 'p_cpi_cpi', [2]);

end;

// Implicit parameters associated to the adl command must be calibrated after the model block.
p_ffr_ffr_lag_1 = 1;
p_ffr_ffr_lag_2 = p_ffr_ffr_lag_1*.5;
p_ffr_ffr_lag_3 = p_ffr_ffr_lag_2*.5;
p_ffr_unrate_lag_1 = 2;
p_ffr_cpi_lag_4 = 3;

// Actual paths for the variables. You migth instead do a stochastic simulation of the model
// to define these paths.
ds1 = dseries(randn(30, 3), 1990Q1, {'ffr', 'unrate', 'cpi'});

// Baseline paths for the variables.
ds0 = dseries(zeros(30, 3), 1990Q1, {'ffr', 'unrate', 'cpi'});

/* Trigger the decomposition in levels of an equation
**
** - First argument (ffr) is the name of the equation to be decomposed.
** - Second argument (ds1) is a dseries object containing the actual paths of the endogenous and exogenous variables.
** - Third argument (ds0) is a dseries object containing the baseline paths of the endogenous and exogenous variables.
**
** If there is no error (missing variables, undefined equations, ...) a figure will pop up with displaying the
** contributions of the variables appearing on the right hand side of equation `ffr`. Note that the overall contribution
** of each variable (at all lags) is reported. If you want to decompose these aggregates, you need to rewrite the
** model as shown in the other mod file in the same folder (contrib2.mod).
*/

plot_contributions ffr ds1 ds0 ;
