conditional_forecast_paths;
var P_obs;
periods  1  2  3:5;
values   1.01 0.98 1;
var Y_obs;
periods  1  2  3:5;
values   3 3.04 2.98;
end;

conditional_forecast(periods=100,parameter_set=posterior_mode,replic=1000, controlled_varexo=(e_y,e_p));

plot_conditional_forecast(periods=100) Y_obs P_obs;

%compare unconditional forecasts
cond_forecast=oo_.conditional_forecast;
if max(abs(cond_forecast.uncond.Mean.P_obs(2:end)-oo_.forecast.Mean.P_obs))>1e-8
    error('Unconditional Forecasts do not match')
end
        
%compare conditional forecasts; histval here sets initval condition for capital different from steady state
initial_condition_states=NaN(size(oo_.dr.ys));
initial_condition_states(strmatch('Y_obs',M_.endo_names,'exact')) = oo_.SmoothedVariables.Y_obs(end)-oo_.Smoother.Trend.Y_obs(end);
initial_condition_states(strmatch('P_obs',M_.endo_names,'exact')) = oo_.SmoothedVariables.P_obs(end)-oo_.Smoother.Trend.P_obs(end);
initial_condition_states(strmatch('junk1',M_.endo_names,'exact')) = oo_.SmoothedVariables.junk1(end);
initial_condition_states(strmatch('junk2',M_.endo_names,'exact')) = oo_.SmoothedVariables.junk2(end)-oo_.Smoother.Trend.junk2(end);
shock_matrix = zeros(options_cond_fcst_.periods ,M_.exo_nbr); %create shock matrix with found controlled shocks
shock_matrix(1:5,strmatch('e_y',M_.exo_names,'exact')) = cond_forecast.controlled_exo_variables.Mean.e_y; %set controlled shocks to their values
shock_matrix(1:5,strmatch('e_p',M_.exo_names,'exact')) = cond_forecast.controlled_exo_variables.Mean.e_p; %set controlled shocks to their values

y_simult = simult_(M_,options_,initial_condition_states,oo_.dr,shock_matrix,1);

if max(abs(y_simult(strmatch('Y_obs',M_.endo_names,'exact'),:)'+(options_.first_obs-1+options_.nobs:options_.first_obs-1+options_.nobs+options_.forecast)'*g_y-cond_forecast.cond.Mean.Y_obs))>1e-8
    error('Conditional Forecasts do not match')
end

if max(abs(y_simult(strmatch('P_obs',M_.endo_names,'exact'),:)'+(options_.first_obs-1+options_.nobs:options_.first_obs-1+options_.nobs+options_.forecast)'*g_p-cond_forecast.cond.Mean.P_obs))>1e-8
    error('Conditional Forecasts do not match')
end

