// --+ options: json=compute +--

/* REMARKS
** -------
**
** Please, firt read the content of contrib1.mod.
**
** I consider the same model as in contrib1.mod, but here we want to obtain the contributions of
** the variables to `ffr` at all lags (individually). To do this we need to add variables for all
** the lags, add corresponding equations (defining the new variables) and add the new variables
** in the dseries objects.
*/

// Declaration of the endogenous variables
var ffr, unrate, cpi;

// Define additional variables
var ffr_lag_1, ffr_lag_2, ffr_lag_3, unrate_lag_1, cpi_lag_4;

// Declaration of the exogenous variables
varexo e_ffr, e_unrate, e_cpi;


// Declare the parameters appearing in the equation to be decomposed.
parameters p_ffr_ffr_lag_1 p_ffr_ffr_lag_2 p_ffr_ffr_lag_3 p_ffr_unrate_lag_1 p_ffr_cpi_lag_4;

// Declaration of the model. Note that you must associate a name to the equations. This is mandatory to
// select the equation for which you need to perform the decomposition.
model;

[name='ffr']
     //ffr = adl(ffr, 'p_ffr_ffr', [1:3]) + adl(unrate, 'p_ffr_unrate', 1) + adl(cpi, 'p_ffr_cpi', [4]);
     ffr = p_ffr_ffr_lag_1*ffr_lag_1 + p_ffr_ffr_lag_2*ffr_lag_2 + p_ffr_ffr_lag_3*ffr_lag_3
           + p_ffr_unrate_lag_1*unrate_lag_1
           + p_ffr_cpi_lag_4*cpi_lag_4;

[name='unrate']
     unrate = adl(unrate, 'p_ffr_unrate', [4 2 5]) + adl(cpi, 'p_unrate_cpi', 6);

[name='cpi']
     cpi = adl(ffr, 'p_cpi_ffr', 2) + adl(cpi, 'p_cpi_cpi', [2]);

// Definitions of the auxiliary variables (we don't need names for these equations).
ffr_lag_1 = ffr(-1);
ffr_lag_2 = ffr(-2);
ffr_lag_3 = ffr(-3);
unrate_lag_1 = unrate(-1);
cpi_lag_4 = cpi(-4);

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

// Create auxiliary variables as dseries objects
verbatim;

  ffr_lag_1 = dseries(ds1.ffr(-1).data, ds1.firstdate, 'ffr_lag_1');
  ffr_lag_2 = dseries(ds1.ffr(-2).data, ds1.firstdate, 'ffr_lag_2');
  ffr_lag_3 = dseries(ds1.ffr(-3).data, ds1.firstdate, 'ffr_lag_3');
  unrate_lag_1 = dseries(ds1.unrate(-1).data, ds1.firstdate, 'unrate_lag_1');
  cpi_lag_4 = dseries(ds1.cpi(-4).data, ds1.firstdate, 'cpi_lag_4');

end;

// Put them in ds1
ds1 = [ds1, ffr_lag_1, ffr_lag_2, ffr_lag_3, unrate_lag_1, cpi_lag_4];

// Baseline paths for the variables.
ds0 = dseries(zeros(30, 8), 1990Q1, {'ffr', 'unrate', 'cpi', 'ffr_lag_1', 'ffr_lag_2', 'ffr_lag_3', 'unrate_lag_1', 'cpi_lag_4'});

/* Trigger the decomposition in levels of an equation
**
** - First argument (ffr) is the name of the equation to be decomposed.
** - Second argument (ds1) is a dseries object containing the actual paths of the endogenous and exogenous variables.
** - Third argument (ds0) is a dseries object containing the baseline paths of the endogenous and exogenous variables.
**
** If there is no error (missing variables, undefined equations, ...) a figure will pop up with displaying the
** contributions of the variables appearing on the right hand side of equation `ffr`.
*/

plot_contributions ffr ds1 ds0 ;
