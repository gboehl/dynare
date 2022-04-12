@#include "../Trend_model_prefilter_common.inc"

addpath('..');
generate_trend_stationary_AR1(M_.fname);

estimation(order=1,datafile='Trend_prefilter_first_obs_AR1_trend_data_with_constant',
    mh_replic=0,mode_compute=4,
    filtered_vars, filter_step_ahead = [1,2,4],        
    first_obs=1000,diffuse_filter,smoother,forecast=100,prefilter=1) P_obs Y_obs junk2;

load('Trend_prefilter_first_obs_AR1_trend_data_with_constant');
@#include "../Trend_load_data_common.inc" 


loaded_par=load('Trend_prefilter_first_obs_orig_params');

if max(abs((M_.params-loaded_par.orig_params([1:4,7:8]))./loaded_par.orig_params([1:4,7:8])))>0.03
    error('Parameter estimates do not match')
end

y_forecast_100_periods=loaded_par.orig_params(strmatch('const_y',loaded_par.param_names,'exact'))+(options_.first_obs+options_.nobs-1+options_.forecast)*loaded_par.orig_params(strmatch('g_y',loaded_par.param_names,'exact'))
p_forecast_100_periods=loaded_par.orig_params(strmatch('const_p',loaded_par.param_names,'exact'))+(options_.first_obs+options_.nobs-1+options_.forecast)*loaded_par.orig_params(strmatch('g_p',loaded_par.param_names,'exact'))

@#include "../Trend_diagnostics_ML_common.inc" 
