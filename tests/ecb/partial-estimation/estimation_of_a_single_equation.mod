// --+ options: json=compute +--

// Declare the endogenous variables.
var ffr, unrate, cpi;

// Declare the exogenous variables.
varexo e_ffr, e_unrate, e_cpi;

// Declare the model.
model(linear);

[name='ffr']
     ffr = adl(ffr, 'p_ffr_ffr', [1:3]) + adl(unrate, 'p_ffr_unrate', 1) + adl(cpi, 'p_ffr_cpi', [4]) + e_ffr;

[name='unrate']
     unrate = adl(unrate, 'p_ffr_unrate', [4 2 5]) + adl(cpi, 'p_unrate_cpi', 6) + e_unrate;

[name='cpi']
     cpi = adl(ffr, 'p_cpi_ffr', 2) + adl(cpi, 'p_cpi_cpi', [2]) + e_cpi;

end;

// Set true values for the parameters (randomly, I write directly in M_.params
// because I am too lazy to write all the parameter names).
M_.params = .1*randn(M_.param_nbr, 1);

shocks;
var e_ffr; stderr .01;
var e_unrate; stderr .01;
var e_cpi; stderr .01;
end;

// Create a dataset by simulating the model.
oo_ = simul_backward_model([], 1000, options_, M_, oo_);

// Put all the simulated data in a dseries object.
ds1 = dseries(transpose(oo_.endo_simul), 1900Q1, cellstr(M_.endo_names));


// Select a subsample for estimation
ds1 = ds1(1951Q1:1990Q1);


// Print true values of the parameters (for comparison with the estimation results).
dyn_table('True values for the parameters (DGP)', {}, {}, cellstr(M_.param_names), ...
    {'Parameter Value'}, 4, M_.params);


/* Run the estimation of equation ffr and display the results. Note that the corresponding
** values in M_.params will be updated.
**
** estimate command
** ----------------
**
** First argument is used to specify the estimation method
** Second argument is used to specify the dataset (an existing dseries object).
** Third argument (and following) is(are) the name(s) of equation(s) to be estimated.
*/

estimate method(ols) data(ds1) ffr;
