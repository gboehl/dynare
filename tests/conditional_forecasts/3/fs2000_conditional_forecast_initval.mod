/* Mod file tests the correctness of the conditional_forecast command when used together with initval by
 * - checking whether the unconditional forecast from the conditional_forecast-command 
 *      coincides with the one from the forecast command
 * - checking whether the conditional forecast coincides with the path of 
 *      capital derived when simulating the model with simult_ and the computed exogenous instruments
 * - initval should not play a role as initial conditions need to be set with histval; with initval
 *      the subsequently computed steady state is used
 */

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

stoch_simul(irf=0);

conditional_forecast_paths;
var gy_obs;
periods  1  2  3:5;
values   0.01 -0.02 0;
var gp_obs;
periods  1  2  3:5;
values   1 1.04 0.98;
end;

%set capital to non-steady state value
initval;
k = 6; 
end;

conditional_forecast(periods=100,parameter_set=calibration,replic=10000, controlled_varexo=(e_m,e_a));

plot_conditional_forecast(periods=100) gy_obs gp_obs k;
forecast(periods=100);

%compare unconditional forecasts
cond_forecast=oo_.conditional_forecast;
if max(abs(cond_forecast.uncond.Mean.k(2:end)-oo_.forecast.Mean.k))>1e-8
    error('Unconditional Forecasts do not match')
end
        
%compare conditional forecasts; initval should not play a role as initial
%conditions need to be set with histval; 
initial_condition_states = oo_.dr.ys; 
shock_matrix = zeros(options_cond_fcst_.periods ,M_.exo_nbr); %create shock matrix with found controlled shocks
shock_matrix(1:5,strmatch('e_a',M_.exo_names,'exact')) = cond_forecast.controlled_exo_variables.Mean.e_a; %set controlled shocks to their values
shock_matrix(1:5,strmatch('e_m',M_.exo_names,'exact')) = cond_forecast.controlled_exo_variables.Mean.e_m; %set controlled shocks to their values

y_simult = simult_(M_,options_,initial_condition_states,oo_.dr,shock_matrix,1);

if max(abs(y_simult(strmatch('k',M_.endo_names,'exact'),:)'-cond_forecast.cond.Mean.k))>1e-8
    error('Unconditional Forecasts do not match')
end
