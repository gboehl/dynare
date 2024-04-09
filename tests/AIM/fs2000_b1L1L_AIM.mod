// This file replicates the estimation of the CIA model from 
// Frank Schorfheide (2000) "Loss function-based evaluation of DSGE models" 
// Journal of  Applied Econometrics, 15, 645-670.
// the data are the ones provided on Schorfheide's web site with the programs.
// http://www.econ.upenn.edu/~schorf/programs/dsgesel.ZIP
// You need to have fsdat.m in the same directory as this file.
// This file replicates: 
// -the posterior mode as computed by Frank's Gauss programs
// -the parameter mean posterior estimates reported in the paper
// -the model probability (harmonic mean) reported in the paper
// This file was tested with dyn_mat_test_0218.zip
// the smooth shocks are probably stil buggy
//
// The equations are taken from J. Nason and T. Cogley (1994) 
// "Testing the implications of long-run neutrality for monetary business
// cycle models" Journal of Applied Econometrics, 9, S37-S70.
// Note that there is an initial minus sign missing in equation (A1), p. S63.
//
// Michel Juillard, February 2004
options_.usePartInfo=0;
var m P c e W R k d n l gy_obs gp_obs Y_obs P_obs y dA P2 c2;
varexo e_a e_m;

parameters alp bet gam mst rho psi del;

alp = 0.33;
bet = 0.99;
gam = 0.003;
mst = 1.011;
rho = 0.7;
psi = 0.787;
del = 0.02;

model ;
dA = exp(gam+e_a);
log(m) = (1-rho)*log(mst) + rho*log(m(-1))+e_m;
-P/(c(+1)*P(+1)*m)+bet*P(+1)*(alp*exp(-alp*(gam+log(e(+1))))*k^(alp-1)*n(+1)^(1-alp)+(1-del)*exp(-(gam+log(e(+1)))))/(c2(+1)*P2(+1)*m(+1))=0;
W = l/n;
-(psi/(1-psi))*(c*P/(1-n))+l/n = 0;
R = P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(-alp)/W;
1/(c*P)-bet*P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)/(m*l*c(+1)*P(+1)) = 0;
c+k = exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)+(1-del)*exp(-(gam+e_a))*k(-1);
P*c = m;
m-1+d = l;
e = exp(e_a);
y = k(-1)^alp*n^(1-alp)*exp(-alp*(gam+e_a));
gy_obs = dA*y/y(-1);
gp_obs = (P/P(-1))*m(-1)/dA;
Y_obs/Y_obs(-1) = gy_obs;
P_obs/P_obs(-1) = gp_obs;
P2 = P(+1);
c2 = c(+1);
end;

steady_state_model;
  dA = exp(gam);
  gst = 1/dA;
  m = mst;
  khst = ( (1-gst*bet*(1-del)) / (alp*gst^alp*bet) )^(1/(alp-1));
  xist = ( ((khst*gst)^alp - (1-gst*(1-del))*khst)/mst )^(-1);
  nust = psi*mst^2/( (1-alp)*(1-psi)*bet*gst^alp*khst^alp );
  n  = xist/(nust+xist);
  P  = xist + nust;
  k  = khst*n;

  l  = psi*mst*n/( (1-psi)*(1-n) );
  c  = mst/P;
  d  = l - mst + 1;
  y  = k^alp*n^(1-alp)*gst^alp;
  R  = mst/bet;
  W  = l/n;
  ist  = y-c;
  q  = 1 - d;

  e = 1;
  
  gp_obs = m/dA;
  gy_obs = dA;
  Y_obs = 1;
  P_obs = 1;
  P2 = P;
  c2 = c;
end;


shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
end;

steady(nocheck);
 
stoch_simul(aim_solver, order=1, irf=0);
 
benchmark = load(['fs2000_b1L1L' filesep 'Output' filesep 'fs2000_b1L1L_results']);
threshold = 1e-8;
 
if max(max(abs(benchmark.oo_.dr.ghx-oo_.dr.ghx))) > threshold
  error('error in ghx');
elseif max(max(abs(benchmark.oo_.dr.ghu-oo_.dr.ghu))) > threshold
  error('error in ghu');
end
