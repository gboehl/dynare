@#include "../Trend_model_prefilter_common.inc"
        
addpath('..');
generate_trend_stationary_AR1(M_.fname);

estimation(order=1,datafile='Trend_prefilter_MC_AR1_trend_data_with_constant',mh_replic=400,mode_compute=4,
    first_obs=1,smoother,prefilter=1,
    mh_nblocks=1,mh_jscale=1e-4,
    filtered_vars, filter_step_ahead = [1,2,4],
    mcmc_jumping_covariance='Trend_prefilter_MC_MCMC_jump_covar_prefilter',forecast=100) P_obs Y_obs junk2;

load('Trend_prefilter_MC_AR1_trend_data_with_constant');
@#include "../Trend_load_data_common.inc" 

loaded_par=load('Trend_prefilter_MC_orig_params_prefilter');

if max(abs((M_.params-loaded_par.orig_params)./loaded_par.orig_params))>0.03
    error('Parameter estimates do not match')
end
loaded_par_full=load('Trend_prefilter_MC_orig_params');
y_forecast_100_periods=loaded_par_full.orig_params(strmatch('const_y',loaded_par_full.param_names,'exact'))+(options_.first_obs+options_.nobs-1+options_.forecast)*loaded_par_full.orig_params(strmatch('g_y',loaded_par_full.param_names,'exact'));
p_forecast_100_periods=loaded_par_full.orig_params(strmatch('const_p',loaded_par_full.param_names,'exact'))+(options_.first_obs+options_.nobs-1+options_.forecast)*loaded_par_full.orig_params(strmatch('g_p',loaded_par_full.param_names,'exact'));

@#include "../Trend_diagnostics_MCMC_common.inc"