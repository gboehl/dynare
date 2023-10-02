@#include "../Trend_exp_model_no_prefilter_common.inc" 

addpath('..');
generate_trend_stationary_AR1(M_.fname);

estimation(order=1,datafile='Trend_loglinear_no_prefilter_MC_Exp_AR1_trend_data_with_constant',mh_replic=400,
    mode_compute=4,silent_optimizer,first_obs=1,loglinear,diffuse_filter,smoother,forecast=100,prefilter=0,
    mcmc_jumping_covariance='Trend_loglinear_no_prefilter_MC_MCMC_jump_covar',
    filtered_vars, filter_step_ahead = [1,2,4],
    mh_nblocks=1,mh_jscale=0.3) P_obs Y_obs junk2;
load('Trend_loglinear_no_prefilter_MC_Exp_AR1_trend_data_with_constant');
@#include "../Trend_load_data_common.inc" 

loaded_par=load('Trend_loglinear_no_prefilter_MC_orig_params');
if max(abs((M_.params-loaded_par.orig_params)./loaded_par.orig_params))>0.03
    error('Parameter estimates do not match')
end
y_forecast_100_periods=loaded_par.orig_params(strmatch('const_y',loaded_par.param_names,'exact'))+(options_.first_obs+options_.nobs-1+options_.forecast)*loaded_par.orig_params(strmatch('g_y',loaded_par.param_names,'exact'));
p_forecast_100_periods=loaded_par.orig_params(strmatch('const_p',loaded_par.param_names,'exact'))+(options_.first_obs+options_.nobs-1+options_.forecast)*loaded_par.orig_params(strmatch('g_p',loaded_par.param_names,'exact'));

@#include "../Trend_diagnostics_MCMC_common.inc"