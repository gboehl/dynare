close all

dynare ols.mod;

global options_
options_.noprint = true;

NSIMS = 1000;

calibrated_values = M_.params;
Sigma_e = M_.Sigma_e;

options_.bnlms.set_dynare_seed_to_default = false;

nparampool = length(M_.params);
BETA = zeros(NSIMS, nparampool);
for i=1:NSIMS
    i
    firstobs = rand(3, length(M_.endo_names));
    M_.params = calibrated_values;
    M_.Sigma_e = Sigma_e; 
    simdata = simul_backward_model(dseries(firstobs, dates('1995Q1'), M_.endo_names), 10000);
    simdata = simdata(simdata.dates(5001:6000));
    names = regexp(simdata.name, 'res\w*');
    idxs = find(cellfun(@isempty, names));
    dyn_ols(simdata{idxs}, {}, {'eq7'});
    BETA(i, :) = M_.params';
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