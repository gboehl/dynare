/* Checks that, for order = 2 and k_order_solver, a model with 2 leads
   and the same model with one lead (using auxiliary vars) give the same result (M-file version).
   The counterfactual is fs2000k2_m.mod */

var m m_1 P P_1 c e W R k d n l gy_obs gp_obs y dA AUXv;
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
log(m) = (1-rho)*log(mst) + rho*log(m_1(-1))+e_m;
-P/(c(+1)*P(+1)*m)+AUXv(+1)=0;
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
gp_obs = (P/P_1(-1))*m_1(-1)/dA;
m_1 = m;
P_1 = P;
AUXv = bet*P*(alp*exp(-alp*(gam+log(e)))*k(-1)^(alp-1)*n^(1-alp)+(1-del)*exp(-(gam+log(e))))/(c(+1)*P(+1)*m);
end;

steady_state_model;
  dA = exp(gam);
  gst = 1/dA;
  m = mst;
  m_1=mst;
  khst = ( (1-gst*bet*(1-del)) / (alp*gst^alp*bet) )^(1/(alp-1));
  xist = ( ((khst*gst)^alp - (1-gst*(1-del))*khst)/mst )^(-1);
  nust = psi*mst^2/( (1-alp)*(1-psi)*bet*gst^alp*khst^alp );
  n  = xist/(nust+xist);
  P  = xist + nust;
  P_1 = P;
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
  AUXv = 1/(c*m);
end;


shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
end;

steady;

stoch_simul(order=2,k_order_solver,irf=0);

if ~exist('fs2000k2_m_results.mat','file');
   error('fs2000k2_m must be run first');
end;

oo1 = load('fs2000k2_m_results','oo_');

dr0 = oo1.oo_.dr;
dr = oo_.dr;

ikr = [2:10 1 13:17];
ikc = [1 3 4 2];
ikc2 = [1 3 4 2 9 11 12 10 13 15 16 14 5 7 8 6];
ikc2u = [1 2 5 6 7 8 3 4];

if max(max(abs(dr0.ghx - dr.ghx(ikr,ikc)))) > 1e-12;
   error('error in ghx');
end;
if max(max(abs(dr0.ghu - dr.ghu(ikr,:)))) > 1e-12;
   error('error in ghu');
end;
if max(max(abs(dr0.ghxx - dr.ghxx(ikr,ikc2)))) > 1e-12;
   error('error in ghxx');
end;
if max(max(abs(dr0.ghuu - dr.ghuu(ikr,:)))) > 1e-12;
   error('error in ghuu');
end;
if max(max(abs(dr0.ghxu - dr.ghxu(ikr,ikc2u)))) > 1e-12;
   error('error in ghxu');
end;
if max(max(abs(dr0.ghs2 - dr.ghs2(ikr,:)))) > 1e-12;
   error('error in ghs2');
end;

