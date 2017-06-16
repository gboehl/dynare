var dx dy;
varobs dx dy;

bvar_density(datafile = bvar_sample, first_obs = 20, bvar_prior_flat,
             bvar_prior_train = 10) 2;

bvar_forecast(forecast = 2, bvar_replic = 1000, nobs = 200) 2;

bvar_irf(2,'Cholesky');
bvar_irf(2,'SquareRoot');