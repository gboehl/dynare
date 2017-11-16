close all

dynare panel_var_diff_NB_simulation_test.mod;

NSIMS = 1000;

calibrated_values = M_.params;
Sigma_e = M_.Sigma_e;

options_.bnlms.set_dynare_seed_to_default = false;

M_endo_names_trim = cellstr(M_.endo_names);
nparampool = length(M_.params);
BETA = zeros(NSIMS, nparampool);
for i=1:NSIMS
    i
    firstobs = rand(3, length(M_endo_names_trim));
    M_.params = calibrated_values;
    M_.Sigma_e = Sigma_e; 
    simdata = simul_backward_model(dseries(firstobs, dates('1995Q1'), M_endo_names_trim), 10000);
    simdata = simdata(simdata.dates(5001:6000));
    pooled_ols(simdata, ...
        {'de','u2'}, ...
        {'*_q_yed_ecm_*_q_yed_L1', ...
        '*_q_yed_ecm_u2_stn_L1', ...
        '*_q_yed_*_g_yer_L1', ...
        '*_q_yed_u2_stn_L1', ...
        '*_g_yer_ecm_*_q_yed_L1', ...
        '*_g_yer_ecm_u2_stn_L1', ...
        '*_g_yer_*_q_yed_L1', ...
        '*_g_yer_*_g_yer_L1', ...
        '*_g_yer_u2_stn_L1', ...
        '*_ehic_*_ehic_L1'});
    BETA(i, :) = M_.params';
    oldsim = simdata;
end

mean(BETA)' - calibrated_values

for i=1:nparampool
    figure
    hold on
    title(strrep(M_.param_names(i,:), '_', '\_'));
    histogram(BETA(:,i),50);
    line([calibrated_values(i) calibrated_values(i)], [0 NSIMS/10], 'LineWidth', 2, 'Color', 'r');
    hold off
end