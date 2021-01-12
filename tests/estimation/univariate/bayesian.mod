// --+ options: json=compute +--
// https://www3.nd.edu/~rwilliam/stats1/OLS-Stata9.pdf

var income educ jobexp race;

varexo e;

parameters b1, b2, b3, b4;

model;

 // The following three lines are required because the variables are declared as endogenous, but the equations
 // have no consequences on the estimation of the last equation in the block.
 educ = educ;
 jobexp = jobexp;
 race = race;

 [name='eqols']
 income = b1 + b2*educ + b3*jobexp + b4*race + e;

end;

// Load data.
ds = dseries('data.csv');

// Estimate last equatioon by OLS (ds is updated with the fitted values).
ds = dyn_ols(ds, {}, {'eqols'});

/*
** BAYESIAN ESTIMATION OF 'eqols'
*/

// Set the number of estimated parameters.
n = length(oo_.ols.eqols.beta);

// We consider a kind of non informative prior by assuming that the inverse of the prior variance for β is zero.
V0 = diag(Inf(n, 1));

// We choose the prior mean for β randomly (around the OLS estimate).
beta0 = oo_.ols.eqols.beta + randn(n, 1); 

// Set prior for the variance of the error
s2priormean = 17;
s2df = 1;

// Run bayesian estimation with normal-gamma prior. Because there is no closed
// form expression for the joint posterior marginal distribution of β, a Gibbs
// sampling algorithm is used (the prior for β and the inverse of σ² are independent).

gibbslength = 1000000; // Set the number of iterations in Gibbs 
burnin = 10000;        // Set the number of iterations to be discarded (try to remove the effects of the initial condition).
steps = 10;            // Do not record all iterations (try to remove the dependence between the draws).

ds = olsgibbs(ds, 'eqols', beta0, V0, s2priormean, s2df, gibbslength, burnin, steps, {'eqols', 'eqols_olsgibbs_fit'});

// Since we use a diffuse prior for β, the posterior mean of β should be close to the OLS estimate.
if max(abs(oo_.ols.eqols.beta-oo_.olsgibbs.eqols.posterior.mean.beta))>.1
   error('Something is wrong in the Gibbs sampling routine (univariate model)')
end

// We can plot the posterior distribution of the parameters:
post = oo_.olsgibbs.eqols.posterior.density.b2;
plot(post(:,1), post(:,2), '-k', 'linewidth', 2);
axis tight
box on
hold on
plot([oo_.olsgibbs.eqols.posterior.mean.beta(2), oo_.olsgibbs.eqols.posterior.mean.beta(2)], [0 max(post(:,2))], '-r')
hold off