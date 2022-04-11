@#include "../Trend_exp_model_calib_prefilter_common.inc"
options_.filter_decomposition=1;

addpath('..');
generate_trend_stationary_AR1(M_.fname);

calib_smoother(datafile='Tr_prefil_f_obs_loglin_cal_smoother_Exp_AR1_trend_data_with_constant',prefilter=1,loglinear,first_obs=1000,
//         filter_decomposition,
        filtered_vars, filter_step_ahead = [1,2,4]) P_obs Y_obs junk2;
load('Tr_prefil_f_obs_loglin_cal_smoother_Exp_AR1_trend_data_with_constant');
@#include "../Trend_load_data_common.inc" 

loaded_par=load('Tr_prefil_f_obs_loglin_cal_smoother_orig_params');
if max(abs((M_.params-loaded_par.orig_params([1:4,7:8]))./loaded_par.orig_params([1:4,7:8])))>0.03
    error('Parameters do not match')
end

@#include "../Trend_diagnostics_calib_common.inc"