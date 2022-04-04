@#include "fs2000_model.inc" 

estimated_params;
alp, 0.356;
gam,  0.0085;
del, 0.01;
stderr e_a, 0.035449;
stderr e_m, 0.008862;
corr e_m, e_a, 0;
stderr gp_obs, 1;
stderr gy_obs, 1;
//corr gp_obs, gy_obs,0;
end;

@#define data_file_name="fsdat_simul_uncorr_ME_missing"

%%default
options_.lik_init=1;
estimation(kalman_algo=0,mode_compute=0,order=1,datafile=@{data_file_name},smoother,filter_decomposition,forecast = 8,filtered_vars,filter_step_ahead=[1,3],irf=20) m P c e W R k d y gy_obs;
fval_algo_0=oo_.likelihood_at_initial_parameters;
%%Multivariate Kalman Filter
options_.lik_init=1;
estimation(kalman_algo=1,mode_compute=0,order=1,datafile=@{data_file_name},smoother,filter_decomposition,forecast = 8,filtered_vars,filter_step_ahead=[1,3],irf=20) m P c e W R k d y gy_obs;
fval_algo_1=oo_.likelihood_at_initial_parameters;
SmoothedMeasurementErrors(:,:,1)=cell2mat(struct2cell(oo_.SmoothedMeasurementErrors));
SmoothedShocks(:,:,1)=cell2mat(struct2cell(oo_.SmoothedShocks));
SmoothedVariables(:,:,1)=cell2mat(struct2cell(oo_.SmoothedVariables));

%%Univariate Kalman Filter
options_.lik_init=1;
estimation(kalman_algo=3,mode_compute=0,order=1,datafile=@{data_file_name},smoother,filter_decomposition,forecast = 8,filtered_vars,filter_step_ahead=[1,3],irf=20) m P c e W R k d y gy_obs;
fval_algo_3=oo_.likelihood_at_initial_parameters;
SmoothedMeasurementErrors(:,:,3)=cell2mat(struct2cell(oo_.SmoothedMeasurementErrors));
SmoothedShocks(:,:,3)=cell2mat(struct2cell(oo_.SmoothedShocks));
SmoothedVariables(:,:,3)=cell2mat(struct2cell(oo_.SmoothedVariables));

%%Diffuse Multivariate Kalman Filter
options_.lik_init=1;
estimation(kalman_algo=2,mode_compute=0,datafile=@{data_file_name},smoother,filter_decomposition,forecast = 8,filtered_vars,filter_step_ahead=[1,3],irf=20) m P c e W R k d y gy_obs;
fval_algo_2=oo_.likelihood_at_initial_parameters;
SmoothedMeasurementErrors(:,:,2)=cell2mat(struct2cell(oo_.SmoothedMeasurementErrors));
SmoothedShocks(:,:,2)=cell2mat(struct2cell(oo_.SmoothedShocks));
SmoothedVariables(:,:,2)=cell2mat(struct2cell(oo_.SmoothedVariables));

%%Diffuse univariate Kalman Filter
options_.lik_init=1;
estimation(kalman_algo=4,mode_compute=0,datafile=@{data_file_name},smoother,filter_decomposition,forecast = 8,filtered_vars,filter_step_ahead=[1,3],irf=20) m P c e W R k d y gy_obs;
fval_algo_4=oo_.likelihood_at_initial_parameters;
SmoothedMeasurementErrors(:,:,4)=cell2mat(struct2cell(oo_.SmoothedMeasurementErrors));
SmoothedShocks(:,:,4)=cell2mat(struct2cell(oo_.SmoothedShocks));
SmoothedVariables(:,:,4)=cell2mat(struct2cell(oo_.SmoothedVariables));


%%Multivariate Fast Kalman Filter
options_.lik_init=1;
estimation(kalman_algo=1,fast_kalman_filter,mode_compute=0,order=1,datafile=@{data_file_name},smoother,filter_decomposition,forecast = 8,filtered_vars,filter_step_ahead=[1,3],irf=20) m P c e W R k d y gy_obs;
fval_algo_5=oo_.likelihood_at_initial_parameters;
SmoothedMeasurementErrors(:,:,5)=cell2mat(struct2cell(oo_.SmoothedMeasurementErrors));
SmoothedShocks(:,:,5)=cell2mat(struct2cell(oo_.SmoothedShocks));
SmoothedVariables(:,:,5)=cell2mat(struct2cell(oo_.SmoothedVariables));

%%Multivariate Fast Kalman Filter
options_.lik_init=1;
estimation(kalman_algo=3,fast_kalman_filter,mode_compute=0,order=1,datafile=@{data_file_name},smoother,filter_decomposition,forecast = 8,filtered_vars,filter_step_ahead=[1,3],irf=20) m P c e W R k d y gy_obs;
fval_algo_6=oo_.likelihood_at_initial_parameters;
SmoothedMeasurementErrors(:,:,6)=cell2mat(struct2cell(oo_.SmoothedMeasurementErrors));
SmoothedShocks(:,:,6)=cell2mat(struct2cell(oo_.SmoothedShocks));
SmoothedVariables(:,:,6)=cell2mat(struct2cell(oo_.SmoothedVariables));


if max(max(abs(SmoothedMeasurementErrors-repmat(SmoothedMeasurementErrors(:,:,1),[1,1,6]))))>1e-8
    error('SmoothedMeasurementErrors do not match')
end

if max(max(abs(SmoothedShocks-repmat(SmoothedShocks(:,:,1),[1,1,6]))))>1e-8
    error('SmoothedShocks do not match')
end

if max(max(abs(SmoothedVariables-repmat(SmoothedVariables(:,:,1),[1,1,6]))))>1e-8
    error('SmoothedVariables do not match')
end

if max(abs([fval_algo_2,fval_algo_3,fval_algo_4,fval_algo_5,fval_algo_6]-fval_algo_1))>1e-8
    error('Likelihoods do not match')
end