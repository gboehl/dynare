// Declaration of the endogenous variables of the DSGE model.
var a g mc mrs n winf pie r rw y;

// Declaration of the exogenous variables of the DSGE model.
varexo e_a e_g e_lam e_ms;

// Declaration of the deep parameters
parameters invsig delta gam rho gampie gamy rhoa rhog bet
	   thetabig omega eps ;

eps=6;
thetabig=2;
bet=0.99;
invsig=2.5;
gampie=1.5;
gamy=0.125;
gam=1;
delta=0.36;
omega=0.54;
rhoa=0.5;
rhog=0.5;
rho=0.5;


model(linear);

	y=y(+1)-(1/invsig)*(r-pie(+1)+g(+1)-g);
	y=a+(1-delta)*n;
	mc=rw+n-y;
	mrs=invsig*y+gam*n-g;
	r=rho*r(-1)+(1-rho)*(gampie*pie+gamy*y)+e_ms;
	rw=rw(-1)+winf-pie;
	a=rhoa*a(-1)+e_a;
	g=rhog*g(-1)+e_g;
	rw=mrs;

	// HYBRID PHILLIPS CURVED USED FOR THE SUMULATIONS:
	//   pie = (omega/(1+omega*bet))*pie(-1)+(bet/(1+omega*bet))*pie(1)+(1-delta)*
      	//   (1-(1-1/thetabig)*bet)*(1-(1-1/thetabig))/((1-1/thetabig)*(1+delta*(eps-1)))/(1+omega*bet)*(mc+e_lam);

	// FORWARD LOOKING PHILLIPS CURVE:
	    pie=bet*pie(+1)+(1-delta)*(1-(1-1/thetabig)*bet)*(1-(1-1/thetabig))/((1-1/thetabig)*(1+delta*(eps-1)))*(mc+e_lam);
end;



// Declaration of the prior beliefs about the deep parameters.
estimated_params;
    stderr e_a, uniform_pdf,,,0,2;
    stderr e_g, uniform_pdf,,,0,2;
    stderr e_ms, uniform_pdf,,,0,2;
    stderr e_lam, uniform_pdf,,,0,2;

    invsig, gamma_pdf, 2.5, 1.76;
    gam, normal_pdf, 1, 0.5;
    rho, uniform_pdf,,,0,1;
    gampie, normal_pdf, 1.5, 0.25;
    gamy, gamma_pdf, 0.125, 0.075;
    rhoa, uniform_pdf,,,0,1;
    rhog, uniform_pdf,,,0,1;
    thetabig, gamma_pdf, 3, 1.42, 1, ;

    // Parameter for the hybrid Phillips curve
    // omega, uniform_pdf,,,0,1;

end;


/*
** Declaration of the observed endogenous variables. Note that they are the variables of the VAR (4 by default) and that we must
** have as many observed variables as exogenous variables.
*/
varobs pie r rw y;

/* REMARK 1.
** The option dsge_var=.8 triggers the estimation of a DSGE-VAR model, with a calibrated dsge prior weight equal to .8.
**
** REMARK 2.
** The option bayesian_irf triggers the computation of the DSGE-VAR and DSGE posterior distribution of the IRFs.
** The Dashed lines are the first, fifth (ie the median) and ninth posterior deciles of the DSGE-VAR's IRFs, the bold dark curve is the
** posterior median of the DSGE's IRfs and the shaded surface covers the space between the first and ninth posterior deciles of the DSGE's IRFs.
*/
estimation(datafile=datarabanal_hybrid,silent_optimizer,first_obs=50,mh_nblocks = 1,nobs=90,dsge_var=.8,optim=('NumgradAlgorithm',3),mode_compute=4,mh_replic=2000,bayesian_irf);
