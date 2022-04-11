@#include "../Trend_exp_model_no_prefilter_common.inc" 

addpath('..');
generate_trend_stationary_AR1(M_.fname);

estimation(order=1,datafile='Trend_loglinear_no_prefilter_first_obs_Exp_AR1_trend_data_with_constant',mh_replic=0,
    mode_compute=4,first_obs=1000,
    filtered_vars, filter_step_ahead = [1,2,4],        
    loglinear,smoother,forecast=100,prefilter=0) P_obs Y_obs junk2;

load('Trend_loglinear_no_prefilter_first_obs_Exp_AR1_trend_data_with_constant');
@#include "../Trend_load_data_common.inc" 

loaded_par=load('Trend_loglinear_no_prefilter_first_obs_orig_params');
if max(abs((M_.params-loaded_par.orig_params)./loaded_par.orig_params))>0.03
    error('Parameter estimates do not match')
end

y_forecast_100_periods=loaded_par.orig_params(strmatch('const_y',M_.param_names,'exact'))+(options_.first_obs+options_.nobs-1+options_.forecast)*loaded_par.orig_params(strmatch('g_y',M_.param_names,'exact'));
p_forecast_100_periods=loaded_par.orig_params(strmatch('const_p',M_.param_names,'exact'))+(options_.first_obs+options_.nobs-1+options_.forecast)*loaded_par.orig_params(strmatch('g_p',M_.param_names,'exact'));

@#include "../Trend_diagnostics_ML_common.inc" 

@#include "../Trend_no_prefilter_conditional_forecast.inc"
