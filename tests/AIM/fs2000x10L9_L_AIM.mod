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

var m P c e W R k d n l gy_obs gp_obs y dA;
varexo e_a e_m;

parameters alp bet gam mst rho psi del;

alp = 0.33;
bet = 0.99;
gam = 0.003;
mst = 1.011;
rho = 0.7;
psi = 0.787;
del = 0.02;

model;
dA = exp(gam+e_a);
log(m) = (1-rho)*log(mst) + rho*log(0.5*m(-1)+0.25*m(-2)+0.13*m(-3)+0.06*m(-4)+0.03*m(-5)+0.015*m(-6)+0.007*m(-7)+0.004*m(-8)+0.003*m(-9)+0.001*m(-10))+e_m;
-P/(((1.3*c(+1)+c(+5)+0.7*c(+9))*(1.3*P(+1)+P(+5)+0.7*P(+9)))*m/9)+bet*((1.3*P(+1)+P(+5)+0.7*P(+9))/3)*(alp*exp(-alp*(gam+log((1.3*e(+1)+e(+5)+0.7*e(+9))/3)))*k^(alp-1)*n(+1)^(1-alp)+(1-del)*exp(-(gam+log(e(+1)))))/(c(+2)*P(+2)*m(+1))=0;
W = l/n;
-(psi/(1-psi))*(c*P/(1-n))+l/n = 0;
R = P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(-alp)/W;
1/(c*P)-bet*P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)/(m*l* (1.3*c(+1)+c(+5)+0.7*c(+9))*(1.3*P(+1)+P(+5)+0.7*P(+9))/9) = 0;
c+k = exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)+(1-del)*exp(-(gam+e_a)*4)*k(-4);
P*c = m;
m-1+d = l;
e = exp(e_a);
y = k(-1)^alp*n^(1-alp)*exp(-alp*(gam+e_a));
gy_obs = dA*y/y(-1);
gp_obs = (P/P(-1))*m(-1)/dA;
end;

initval;
k = 6;
m = mst;
P = 2.25;
c = 0.45;
e = 1;
W = 4;
R = 1.02;
d = 0.85;
n = 0.19;
l = 0.86;
y = 0.6;
gy_obs = exp(gam);
gp_obs = exp(-gam); 
dA = exp(gam);
end;

shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
end;

steady;

check;

stoch_simul(aim_solver, order=1,irf=0);
 
benchmark = load(['fs2000x10L9_L' filesep 'Output' filesep 'fs2000x10L9_L_results']);
threshold = 1e-8;
 
if max(max(abs(benchmark.oo_.dr.ghx-oo_.dr.ghx))) > threshold
  error('error in ghx');
elseif max(max(abs(benchmark.oo_.dr.ghu-oo_.dr.ghu))) > threshold
  error('error in ghu');
end
