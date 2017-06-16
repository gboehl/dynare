//test selected_variables_only option by comparing results with option to selection of estimation without the option
//allows testing feature, because there is a static variable

var Y_obs P_obs junk1 junk_backwards junk2 junk_static ;
varexo e_y e_p eps_junk;
parameters rho_y rho_p  g_y g_p const_y const_p sigma_y sigma_p;

rho_y=0.5;
rho_p=0.5;
g_y=0.0001;
g_p=-0.0001;
const_y=2;
const_p=2;
sigma_y=0.001;
sigma_p=0.001;


model;
Y_obs = (1-rho_y)*const_y + rho_y*Y_obs(-1)+sigma_y*e_y;
P_obs = (1-rho_p)*const_p + rho_p*P_obs(-1)+sigma_p*e_p;
junk1 = 0.9*junk1(+1);
junk2 = 0.9*junk2(-1)+eps_junk;
junk_backwards=0.9*junk_backwards(-1);
junk_static=junk_backwards;
end;

steady_state_model;
Y_obs = const_y;
P_obs = const_p;
junk1=0;
junk2=0;
junk_backwards=0;
end;

shocks;
var e_p; stderr 1;
var e_y; stderr 1;
var eps_junk; stderr 1;
end;

steady;

estimated_params;
const_y, normal_pdf, 2, 1;
const_p, normal_pdf, 2, 1;
g_y, normal_pdf, 0.0001, 1;
g_p, normal_pdf, -0.0001, 1;
rho_y, normal_pdf, 0.5, 1;
rho_p, normal_pdf, 0.5, 1;
sigma_y, inv_gamma_pdf, 0.001, inf;
sigma_p, inv_gamma_pdf, 0.001, inf;
end;

options_.plot_priors=0;
varobs P_obs Y_obs junk2;

observation_trends;
Y_obs (g_y);
P_obs (g_p);
end;

estimated_params_init(use_calibration);
end;

addpath('..');
generate_trend_stationary_AR1;

estimation(order=1,datafile='AR1_trend_data_with_constant',mh_replic=0,
        mode_compute=4,first_obs=1,nobs=1000,
        filtered_vars, filter_step_ahead = [1,2,4],        
        diffuse_filter,filter_covariance,smoother,forecast=0,prefilter=0,filter_decomposition) P_obs Y_obs junk2;

%Test selected_variables_only option
oo_all_variables=oo_;
%reset oo_
oo_.Smoother  = [];
oo_.FilteredVariablesKStepAhead  = [];
oo_.FilteredVariablesKStepAheadVariances  = [];
oo_.SmoothedVariables  = [];
oo_.FilteredVariables  = [];
oo_.UpdatedVariables  = [];
oo_.SmoothedShocks  = [];

set_dynare_seed('default');
estimation(order=1,datafile='AR1_trend_data_with_constant',mh_replic=0,
        mode_compute=4,first_obs=1,nobs=1000,
        filtered_vars, filter_step_ahead = [1,2,4],        
        diffuse_filter,smoother,forecast=0,filter_covariance,prefilter=0,filter_decomposition,selected_variables_only) P_obs Y_obs junk2;

% do checks
        
if max(abs(oo_.SmoothedVariables.Y_obs-oo_all_variables.SmoothedVariables.Y_obs))>1e-8 ||...
    max(abs(oo_.SmoothedVariables.P_obs-oo_all_variables.SmoothedVariables.P_obs))>1e-8 || ...
    max(abs(oo_.SmoothedVariables.junk2-oo_all_variables.SmoothedVariables.junk2))>1e-8
    error('Smoothed Variables are wrong')
end

if max(abs(oo_.UpdatedVariables.Y_obs-oo_all_variables.UpdatedVariables.Y_obs))>1e-8 ||...
    max(abs(oo_.UpdatedVariables.P_obs-oo_all_variables.UpdatedVariables.P_obs))>1e-8 || ...
    max(abs(oo_.UpdatedVariables.junk2-oo_all_variables.UpdatedVariables.junk2))>1e-8
    error('Updated Variables are wrong')
end

if mean(abs(oo_.FilteredVariables.Y_obs-oo_all_variables.FilteredVariables.Y_obs))>1e-8 ||...
    mean(abs(oo_.FilteredVariables.P_obs-oo_all_variables.FilteredVariables.P_obs))>1e-8 
    error('Smoothed Variables are wrong')
end

   
Y_pos=strmatch('Y_obs',M_.endo_names,'exact');
P_pos=strmatch('P_obs',M_.endo_names,'exact');
junk2_pos=strmatch('junk2',M_.endo_names,'exact');

[junk_arg,order_index]=sort([Y_pos,P_pos,junk2_pos]);

   
if max(max(max(abs(oo_.FilteredVariablesKStepAhead-oo_all_variables.FilteredVariablesKStepAhead(:,[Y_pos;P_pos;junk2_pos],:)))))>1e-8
    error('FilteredVariablesKStepAhead is wrong')
end

if options_.filter_covariance && max(max(max(max(abs(oo_.FilteredVariablesKStepAheadVariances-oo_all_variables.FilteredVariablesKStepAheadVariances(:,[Y_pos;P_pos;junk2_pos],[Y_pos;P_pos;junk2_pos],:))))))>1e-8
    error('FilteredVariablesKStepAheadVariances is wrong')
end

if max(max(max(max(abs(oo_.FilteredVariablesShockDecomposition-oo_all_variables.FilteredVariablesShockDecomposition(:,[Y_pos;P_pos;junk2_pos],:,:))))))>1e-8
    error('FilteredVariablesShockDecomposition is wrong')
end
        
if max(abs(oo_.SmoothedShocks.e_y-oo_all_variables.SmoothedShocks.e_y))>1e-8 || max(abs(oo_.SmoothedShocks.e_p-oo_all_variables.SmoothedShocks.e_p))>1e-8
    error('Smoothed Shocks are wrong')
end