/* 
Compare simulation results with pruning at second order
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

T = 10;
ex=randn(T,M_.exo_nbr);

% 2nd order
stoch_simul(order=2, nograph, irf=0);

% simult_.m: implements KKSS
options_.pruning = true;
Y2_simult = simult_(M_,options_,oo_.dr.ys,oo_.dr,ex,options_.order); 

% local_state_space_iteration_2 mex: implements KKSS
constant = oo_.dr.ys(oo_.dr.order_var)+0.5*oo_.dr.ghs2;
ss = oo_.dr.ys(oo_.dr.order_var);
Y1_local = zeros(M_.endo_nbr,T+1);
Y1_local(:,1) = oo_.dr.ys;
Y2_local = zeros(M_.endo_nbr,T+1);
Y2_local(:,1) = oo_.dr.ys;
state_var = oo_.dr.order_var(M_.nstatic+1:M_.nstatic+M_.nspred);
for t=2:T+1
   u = ex(t-1,:)';
   yhat1 = Y1_local(state_var,t-1)-oo_.dr.ys(state_var);
   yhat2 = Y2_local(state_var,t-1)-oo_.dr.ys(state_var);
   [Y2, Y1] = local_state_space_iteration_2(yhat2, u, oo_.dr.ghx, oo_.dr.ghu, constant, oo_.dr.ghxx, oo_.dr.ghuu, oo_.dr.ghxu, yhat1, ss, options_.threads.local_state_space_iteration_2);
   Y1_local(oo_.dr.order_var,t) = Y1;
   Y2_local(oo_.dr.order_var,t) = Y2;
end

if (max(abs(Y2_local(:) - Y2_simult(:)))>1e-12)
   error('2nd-order output of simult_ and local_state_space_iteration_2 output are inconsistent.')
end

% pruned_state_space_system.m: implements Andreasen et al.
pss = pruned_state_space_system(M_, options_, oo_.dr, [], 0, false, false);
Y2_an = zeros(M_.endo_nbr,T+1);
Y2_an(:,1) = oo_.dr.ys;
% z    = [xf;xs;kron(xf,xf)]
% inov = [u;kron(u,u)-E_uu(:);kron(xf,u)]
z = zeros(2*M_.nspred+M_.nspred^2, 1);
z(1:2*M_.nspred, 1) = repmat(Y2_an(state_var,1)-oo_.dr.ys(state_var), 2, 1);
z(2*M_.nspred+1:end, 1) = kron(z(1:M_.nspred, 1), z(M_.nspred+1:2*M_.nspred, 1));
inov = zeros(M_.exo_nbr+M_.exo_nbr^2+M_.nspred*M_.exo_nbr, 1);
for t=2:T+1
   u = ex(t-1,:)';
   inov(1:M_.exo_nbr, 1) = u;
   inov(M_.exo_nbr+1:M_.exo_nbr+M_.exo_nbr^2, 1) = kron(u, u) - M_.Sigma_e(:);
   inov(M_.exo_nbr+M_.exo_nbr^2+1:end, 1) = kron(z(1:M_.nspred, 1), u);
   Y = oo_.dr.ys(oo_.dr.order_var) + pss.d + pss.C*z + pss.D*inov;
   x1 = z(1:M_.nspred,1);
   x2 = z(M_.nspred+1:2*M_.nspred,1);
   z = pss.c + pss.A*z + pss.B*inov;
   Y2_an(oo_.dr.order_var,t) = Y;
end

if (max(abs(Y2_an(:) - Y2_local(:)))>1e-12)
   error('2nd-order output of pruned_state_space_system and local_state_space_iteration_2 output are inconsistent.')
end

% 3rd order
stoch_simul(order=3, nograph, irf=0);
% simult_.m
Y3_simult = simult_(M_,options_,oo_.dr.ys,oo_.dr,ex,options_.order); 
% pruned_state_space_system.m
pss = pruned_state_space_system(M_, options_, oo_.dr, [], 0, false, false);
Y3_an = zeros(M_.endo_nbr,T+1);
Y3_an(:,1) = oo_.dr.ys;
% z    = [xf; xs; kron(xf,xf); xrd; kron(xf,xs); kron(kron(xf,xf),xf)]
% inov = [u; kron(u,u)-E_uu(:); kron(xf,u); kron(xs,u); kron(kron(xf,xf),u); kron(kron(xf,u),u); kron(kron(u,u),u))]
x_nbr    = M_.nspred;
u_nbr    = M_.exo_nbr;
z_nbr    = x_nbr + x_nbr + x_nbr^2 + x_nbr + x_nbr^2 + x_nbr^3;
inov_nbr = u_nbr + u_nbr^2 + x_nbr*u_nbr + x_nbr*u_nbr + x_nbr^2*u_nbr + x_nbr*u_nbr^2 + u_nbr^3;
z = zeros(z_nbr, 1);
inov = zeros(inov_nbr, 1);
for t=2:T+1
   u = ex(t-1,:)';
   xf = z(1:M_.nspred,1);
   xs = z(M_.nspred+1:2*M_.nspred,1);
   inov(1:u_nbr, 1) = u;
   ub = u_nbr;
   inov(ub+1:ub+u_nbr^2, 1) = kron(u, u) - M_.Sigma_e(:);
   ub = ub+u_nbr^2;
   inov(ub+1:ub+x_nbr*u_nbr, 1) = kron(xf, u);
   ub = ub+x_nbr*u_nbr;
   inov(ub+1:ub+x_nbr*u_nbr, 1) = kron(xs, u);
   ub = ub+x_nbr*u_nbr;
   inov(ub+1:ub+x_nbr^2*u_nbr) = kron(kron(xf,xf),u);
   ub = ub+x_nbr^2*u_nbr;
   inov(ub+1:ub+x_nbr*u_nbr^2) = kron(kron(xf,u),u);
   ub = ub+x_nbr*u_nbr^2;
   inov(ub+1:end) = kron(kron(u,u),u);
   Y = oo_.dr.ys(oo_.dr.order_var) + pss.d + pss.C*z + pss.D*inov;
   z = pss.c + pss.A*z + pss.B*inov;
   Y3_an(oo_.dr.order_var,t) = Y;
end

if (max(abs(Y3_an(:) - Y3_simult(:)))>1e-12)
   error('3rd-order outputs of pruned_state_space_system and simult_ are inconsistent.')
end

% local_state_space_iteration_3 mex
Y3_local = zeros(M_.endo_nbr,T+1);
Y3_local(:,1) = oo_.dr.ys;
Ylat_local = zeros(3*M_.endo_nbr,T+1);
Ylat_local(:,1) = repmat(oo_.dr.ys,3,1);
for t=2:T+1
   u = ex(t-1,:)';
   yhat1 = Ylat_local(state_var,t-1)-oo_.dr.ys(state_var);
   yhat2 = Ylat_local(M_.endo_nbr+state_var,t-1)-oo_.dr.ys(state_var);
   yhat3 = Ylat_local(2*M_.endo_nbr+state_var,t-1)-oo_.dr.ys(state_var);
   ylat = [yhat1;yhat2;yhat3];
   [Y3,Y] = local_state_space_iteration_3(ylat, u, oo_.dr.ghx, oo_.dr.ghu,  oo_.dr.ghxx, oo_.dr.ghuu, oo_.dr.ghxu, oo_.dr.ghs2, oo_.dr.ghxxx, oo_.dr.ghuuu, oo_.dr.ghxxu, oo_.dr.ghxuu, oo_.dr.ghxss, oo_.dr.ghuss, ss, options_.threads.local_state_space_iteration_3, true);
   Ylat_local(oo_.dr.order_var,t) = Y(1:M_.endo_nbr);
   Ylat_local(M_.endo_nbr+oo_.dr.order_var,t) = Y(M_.endo_nbr+1:2*M_.endo_nbr);
   Ylat_local(2*M_.endo_nbr+oo_.dr.order_var,t) = Y(2*M_.endo_nbr+1:end);
   Y3_local(oo_.dr.order_var,t) = Y3;
end

Y3_local_1 = Ylat_local(1:M_.endo_nbr,:);
Y3_local_2 = Y3_local_1 + Ylat_local(M_.endo_nbr+1:2*M_.endo_nbr,:) - oo_.dr.ys;

if (max(abs(Y3_local_1(:) - Y1_local(:)))>1e-12)
   error('1st-order outputs of local_state_space_iteration_3 and local_state_space_iteration_2 are inconsistent.')
end

if (max(abs(Y3_local_2(:) - Y2_local(:)))>4e-12)
   error('2nd-order outputs of local_state_space_iteration_3 and local_state_space_iteration_2 are inconsistent.')
end

if (max(abs(Y3_local(:) - Y3_simult(:)))>1e-12)
   error('3rd-order output of simult_ and local_state_space_iteration_3 output are inconsistent.')
end
