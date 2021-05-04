// See fs2000.mod in the examples/ directory for details on the model

var m P c e W R k d n l gy_obs gp_obs y dA;
varexo e_a e_m e_b;

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
-P/(c(+1)*P(+1)*m)+bet*exp(e_b)*P(+1)*(alp*exp(-alp*(gam+log(e(+1))))*k^(alp-1)*n(+1)^(1-alp)+(1-del)*exp(-(gam+log(e(+1)))))/(c(+2)*P(+2)*m(+1))=0;
W = l/n;
-(psi/(1-psi))*(c*P/(1-n))+l/n = 0;
R = P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(-alp)/W;
1/(c*P)-bet*exp(e_b)*P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)/(m*l*c(+1)*P(+1)) = 0;
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

write_latex_dynamic_model;

shocks;
var e_a; stderr 0.014;
var e_b; stderr 0.1;
var e_m; stderr 0.005;
end;

steady;

check;


varobs gp_obs gy_obs;

calib_smoother(datafile=fsdat_simul, filtered_vars, filter_step_ahead = [3:4],filter_covariance,smoothed_state_uncertainty) m P c e W R k d n l y dA;
oo0=oo_;

calib_smoother(datafile=fsdat_simul, filtered_vars, filter_step_ahead = [3:4],filter_covariance,smoothed_state_uncertainty,smoother_redux) m P c e W R k d n l y dA;
oo1=oo_;

calib_smoother(datafile=fsdat_simul, filtered_vars, filter_step_ahead = [3:4],filter_covariance,smoothed_state_uncertainty,kalman_algo=2,smoother_redux) m P c e W R k d n l y dA;
oo2=oo_;


for k=1:M_.exo_nbr
    mserr(k)=max(abs(oo0.SmoothedShocks.(M_.exo_names{k})-oo1.SmoothedShocks.(M_.exo_names{k})));
end
if max(mserr)>1.e-12
    error('smoother_redux with kalman_algo=1 does not replicate original smoother for shocks!')
end
        
vlist = M_.endo_names(oo_.dr.order_var(oo_.dr.restrict_var_list));
for k=1:length(vlist)
    merr(k)=max(abs(oo0.SmoothedVariables.(vlist{k})-oo1.SmoothedVariables.(vlist{k})));
    merrU(k)=max(abs(oo0.UpdatedVariables.(vlist{k})-oo1.UpdatedVariables.(vlist{k})));
    merrF(k)=max(abs(oo0.FilteredVariables.(vlist{k})-oo1.FilteredVariables.(vlist{k})));
    
end
if max(merr)>1.e-12
    error('smoother_redux with kalman_algo=1 does not replicate original smoothed restricted var list!')
end
if max(merrU)>1.e-12
    error('smoother_redux with kalman_algo=1 does not replicate original updated restricted var list!')
end
if max(merrF)>1.e-12
    error('smoother_redux with kalman_algo=1 does not replicate original filtered restricted var list!')
end
        
vlist1 = M_.endo_names(~ismember(M_.endo_names,vlist));
for k=1:length(vlist1)
    merr1(k)=max(abs(oo0.SmoothedVariables.(vlist1{k})-oo1.SmoothedVariables.(vlist1{k})));
    merr1U(k)=max(abs(oo0.UpdatedVariables.(vlist1{k})-oo1.UpdatedVariables.(vlist1{k})));
    merr1F(k)=max(abs(oo0.FilteredVariables.(vlist1{k})-oo1.FilteredVariables.(vlist1{k})));
end
if max(merr1)>1.e-12
    for k=1:length(vlist1)
        merr2(k)=max(abs(oo0.SmoothedVariables.(vlist1{k})(2:end)-oo1.SmoothedVariables.(vlist1{k})(2:end)));
    end
    if max(merr2)>1.e-12
        error('smoother_redux with kalman_algo=1 does not replicate original smoothed static variables!')
    end
end
if max(merr1U)>1.e-12
    for k=1:length(vlist1)
        merr2U(k)=max(abs(oo0.UpdatedVariables.(vlist1{k})(2:end)-oo1.UpdatedVariables.(vlist1{k})(2:end)));
    end
    if max(merr2U)>1.e-12
        error('smoother_redux with kalman_algo=1 does not replicate original updated static variables!')
    end
end
if max(merr1F)>1.e-12
    for k=1:length(vlist1)
        merr2F(k)=max(abs(oo0.FilteredVariables.(vlist1{k})(2:end)-oo1.FilteredVariables.(vlist1{k})(2:end)));
    end
    if max(merr2F)>1.e-12
        error('smoother_redux with kalman_algo=1 does not replicate original filtered static variables!')
    end
end
merrK = max(max(max(abs(oo0.FilteredVariablesKStepAhead-oo1.FilteredVariablesKStepAhead))));
if max(merrK)>1.e-12
    error('smoother_redux with kalman_algo=1 does not replicate original k-step ahead forecasts!')
end
verrK = max(max(max(max(abs(oo0.FilteredVariablesKStepAheadVariances(:,[1:14 16],[1:14 16],:)-oo1.FilteredVariablesKStepAheadVariances(:,[1:14 16],[1:14 16],:))))));         
if verrK>1.e-12
    error('smoother_redux with kalman_algo=1 does not replicate original k-step ahead forecast variances!')
end
verr=max(max(max(abs(oo0.Smoother.Variance([1:14 16],[1:14 16],:)-oo1.Smoother.Variance([1:14 16],[1:14 16],:)))));
if verr>1.e-12
    error('smoother_redux with kalman_algo=1 does not replicate original filter covariance!')
end
verrS=max(max(max(abs(oo0.Smoother.State_uncertainty([1:14 16],[1:14 16],:)-oo1.Smoother.State_uncertainty([1:14 16],[1:14 16],:)))));
if verrS>1.e-12
    error('smoother_redux with kalman_algo=1 does not replicate original state covariance!')
end

// now I check kalman_algo=2
for k=1:M_.exo_nbr
    mserr(k)=max(abs(oo0.SmoothedShocks.(M_.exo_names{k})-oo2.SmoothedShocks.(M_.exo_names{k})));
end
if max(mserr)>1.e-12
    error('smoother_redux with kalman_algo=2 does not replicate original smoother for shocks!')
end
        
vlist = M_.endo_names(oo_.dr.order_var(oo_.dr.restrict_var_list));
for k=1:length(vlist)
    merr(k)=max(abs(oo0.SmoothedVariables.(vlist{k})-oo2.SmoothedVariables.(vlist{k})));
    merrU(k)=max(abs(oo0.UpdatedVariables.(vlist{k})-oo2.UpdatedVariables.(vlist{k})));
    merrF(k)=max(abs(oo0.FilteredVariables.(vlist{k})-oo2.FilteredVariables.(vlist{k})));
    
end
if max(merr)>1.e-12
    error('smoother_redux with kalman_algo=2 does not replicate original smoothed restricted var list!')
end
if max(merrU)>1.e-12
    error('smoother_redux with kalman_algo=2 does not replicate original updated restricted var list!')
end
if max(merrF)>1.e-12
    error('smoother_redux with kalman_algo=2 does not replicate original filtered restricted var list!')
end
        
vlist1 = M_.endo_names(~ismember(M_.endo_names,vlist));
for k=1:length(vlist1)
    merr1(k)=max(abs(oo0.SmoothedVariables.(vlist1{k})-oo2.SmoothedVariables.(vlist1{k})));
    merr1U(k)=max(abs(oo0.UpdatedVariables.(vlist1{k})-oo2.UpdatedVariables.(vlist1{k})));
    merr1F(k)=max(abs(oo0.FilteredVariables.(vlist1{k})-oo2.FilteredVariables.(vlist1{k})));
end
if max(merr1)>1.e-12
    for k=1:length(vlist1)
        merr2(k)=max(abs(oo0.SmoothedVariables.(vlist1{k})(2:end)-oo2.SmoothedVariables.(vlist1{k})(2:end)));
    end
    if max(merr2)>1.e-12
        error('smoother_redux with kalman_algo=2 does not replicate original smoothed static variables!')
    end
end
if max(merr1U)>1.e-12
    for k=1:length(vlist1)
        merr2U(k)=max(abs(oo0.UpdatedVariables.(vlist1{k})(2:end)-oo2.UpdatedVariables.(vlist1{k})(2:end)));
    end
    if max(merr2U)>1.e-12
        error('smoother_redux with kalman_algo=2 does not replicate original updated static variables!')
    end
end
if max(merr1F)>1.e-12
    for k=1:length(vlist1)
        merr2F(k)=max(abs(oo0.FilteredVariables.(vlist1{k})(2:end)-oo2.FilteredVariables.(vlist1{k})(2:end)));
    end
    if max(merr2F)>1.e-12
        error('smoother_redux with kalman_algo=2 does not replicate original filtered static variables!')
    end
end
merrK = max(max(max(abs(oo0.FilteredVariablesKStepAhead-oo2.FilteredVariablesKStepAhead))));
if max(merrK)>1.e-12
    error('smoother_redux with kalman_algo=2 does not replicate original k-step ahead forecasts!')
end
verrK = max(max(max(max(abs(oo0.FilteredVariablesKStepAheadVariances(:,[1:14 16],[1:14 16],:)-oo2.FilteredVariablesKStepAheadVariances(:,[1:14 16],[1:14 16],:))))));         
if verrK>1.e-12
    error('smoother_redux with kalman_algo=2 does not replicate original k-step ahead forecast variances!')
end
verr=max(max(max(abs(oo0.Smoother.Variance([1:14 16],[1:14 16],:)-oo2.Smoother.Variance([1:14 16],[1:14 16],:)))));
if verr>1.e-12
    error('smoother_redux with kalman_algo=2 does not replicate original filter covariance!')
end
verrS=max(max(max(abs(oo0.Smoother.State_uncertainty([1:14 16],[1:14 16],:)-oo2.Smoother.State_uncertainty([1:14 16],[1:14 16],:)))));
if verrS>1.e-12
    error('smoother_redux with kalman_algo=2 does not replicate original state covariance!')
end       
