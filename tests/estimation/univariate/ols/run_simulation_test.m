close all

dynare ols.mod;

global M_ options_ oo_


options_.noprint = true;

number_of_simulations = 1000;

calibrated_values = M_.params;
Sigma_e = M_.Sigma_e;

options_.bnlms.set_dynare_seed_to_default = false;

Beta = zeros(number_of_simulations, M_.param_nbr);

for i=1:number_of_simulations
    % Set initial conditions randomly
    firstobs = rand(3, length(M_.endo_names));
    % Set parameters to calibrated values (because after the
    % estimation parameters in equation 7 are updated with OLS
    % estimator).
    M_.params = calibrated_values;
    M_.Sigma_e = Sigma_e;
    % Simulate the model.
    simdata = simul_backward_model(dseries(firstobs, dates('1995Q1'), M_.endo_names), 10000);
    % Select a subsample.
    simdata = simdata(simdata.dates(5001:6000));
    % Perform the estimation of equation 7.
    names = regexp(simdata.name, 'res\w*');
    idxs = find(cellfun(@isempty, names));
    dyn_ols(simdata{idxs}, {}, {'eq7'});
    % Store the estimation results in Beta
    Beta(i, :) = M_.params';
end

pid = oo_.ols.eq7.param_idxs;

for i=1:length(pid)
    figure(i)
    hold on
    title(strrep(M_.param_names(pid(i),:), '_', '\_'));
    bandwidth = mh_optimal_bandwidth(Beta(:,pid(i)), length(Beta(:,pid(i))), -1, 'gaussian');
    [abscissa, f] = kernel_density_estimate(Beta(:,pid(i)), 256, length(Beta(:,pid(i))), bandwidth, 'gaussian');
    plot(abscissa, f, '-k', 'linewidth', 2);
    line([calibrated_values(pid(i)) calibrated_values(pid(i))], [0 max(f)*1.05], 'LineWidth', 2, 'Color', 'r');
    axis tight
    box on
    hold off
end