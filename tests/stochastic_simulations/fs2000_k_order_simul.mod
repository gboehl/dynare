/* 
Compare simulation results from mex-file and Matlab without pruning at second order
 */

var m P c e W R k d n l gy_obs gp_obs y dA ;
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
log(m) = (1-rho)*log(mst) + rho*log(m(-1))+e_m;
-P/(c(+1)*P(+1)*m)+bet*P(+1)*(alp*exp(-alp*(gam+log(e(+1))))*k^(alp-1)*n(+1)^(1-alp)+(1-del)*exp(-(gam+log(e(+1)))))/(c(+2)*P(+2)*m(+1))=0;
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
end;

shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
end;

steady;

stoch_simul(order=2,k_order_solver,nograph);

ex=randn(5,M_.exo_nbr);

Y2_matlab = simult_(M_,options_,oo_.dr.ys,oo_.dr,ex,options_.order);
Y2_mex = k_order_simul(options_.order,M_.nstatic,M_.npred,M_.nboth,M_.nfwrd,M_.exo_nbr,   ...
                       oo_.dr.ys(oo_.dr.order_var),ex',oo_.dr.ys(oo_.dr.order_var),oo_.dr,...
                       options_.pruning);
Y2_mex(oo_.dr.order_var,:) = Y2_mex;

if max(abs(Y2_matlab(:)-Y2_mex(:)))>1e-8;
   error('2nd order: Matlab and mex simulation routines do not return similar results')
end

Y2_mexiter = NaN(M_.endo_nbr,size(ex,1)+1);
Y2_mexiter(:,1) = oo_.dr.ys;
for it = 1:size(ex,1)
    Y2_temp = k_order_simul(options_.order,M_.nstatic,M_.npred,M_.nboth,M_.nfwrd, ...
                            M_.exo_nbr,Y2_mexiter(oo_.dr.order_var,it),ex(it,:)', ...
                            oo_.dr.ys(oo_.dr.order_var),oo_.dr,options_.pruning);
    Y2_temp(oo_.dr.order_var,:) = Y2_temp;
    Y2_mexiter(:,it+1) = Y2_temp(:,2);    
end

if max((abs(Y2_mex(:) - Y2_mexiter(:))))>1e-8;
   error('2nd order: sequential call does not return similar results')
end

stoch_simul(order=3,k_order_solver,nograph,irf=0);
Y3_mex = k_order_simul(options_.order,M_.nstatic,M_.npred,M_.nboth,M_.nfwrd,M_.exo_nbr,   ...
                       oo_.dr.ys(oo_.dr.order_var),ex',oo_.dr.ys(oo_.dr.order_var),oo_.dr,...
                       options_.pruning);
Y3_mex(oo_.dr.order_var,:) = Y3_mex;
Y3_mexiter = NaN(M_.endo_nbr,size(ex,1)+1);
Y3_mexiter(:,1) = oo_.dr.ys;
for it = 1:size(ex,1)
    Y3_temp = k_order_simul(options_.order,M_.nstatic,M_.npred,M_.nboth,M_.nfwrd, ...
                            M_.exo_nbr,Y3_mexiter(oo_.dr.order_var,it),ex(it,:)', ...
                            oo_.dr.ys(oo_.dr.order_var),oo_.dr,options_.pruning);
    Y3_temp(oo_.dr.order_var,:) = Y3_temp;
    Y3_mexiter(:,it+1) = Y3_temp(:,2);    
end

if max((abs(Y3_mex(:) - Y3_mexiter(:))))>1e-8;
   error('3rd order: sequential call does not return similar results')
end

stoch_simul(order=2,k_order_solver,pruning,nograph,irf=0);

Y2_matlab = simult_(M_,options_,oo_.dr.ys,oo_.dr,ex,options_.order); 
Y2_mex = k_order_simul(options_.order,M_.nstatic,M_.npred,M_.nboth,M_.nfwrd,M_.exo_nbr,   ...
                       oo_.dr.ys(oo_.dr.order_var),ex',oo_.dr.ys(oo_.dr.order_var),oo_.dr,...
                       options_.pruning);
Y2_mex(oo_.dr.order_var,:) = Y2_mex;

if max(abs(Y2_matlab(:)-Y2_mex(:)))>1e-8;
   error('2nd order with pruning: Matlab and mex simulation routines do not return similar results')
end

stoch_simul(order=3,k_order_solver,pruning,nograph,irf=0);

Y3_matlab = simult_(M_,options_,oo_.dr.ys,oo_.dr,ex,options_.order); 
Y3_mex = k_order_simul(options_.order,M_.nstatic,M_.npred,M_.nboth,M_.nfwrd,M_.exo_nbr,   ...
                       oo_.dr.ys(oo_.dr.order_var),ex',oo_.dr.ys(oo_.dr.order_var),oo_.dr,...
                       options_.pruning);
Y3_mex(oo_.dr.order_var,:) = Y3_mex;

if max(abs(Y3_matlab(:)-Y3_mex(:)))>1e-8;
   error('3rd order with pruning: Matlab and mex simulation routines do not return similar results')
end

