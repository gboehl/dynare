clearvars
clearvars -global
close all

!rm nkm_saved_data.mat
!rm my_var_est.mat
!rm nkm_var_saved_data.mat

%% Call Dynare without VAR
dynare example1.mod
save('nkm_saved_data.mat')

%% Estimate VAR using simulated data
disp('VAR Estimation');


% MLS estimate of mu and B (autoregressive_matrices)
% Y = mu + B*Z
% from New Introduction to Multiple Time Series Analysis
Y = oo_.endo_simul(1:2, 2:end);
Z = [ ...
    ones(1, size(Y,2)); ...
    oo_.endo_simul(1:2, 1:end-1); ...
    ];
%B = Y*Z'*inv(Z*Z');
B = Y*Z'/(Z*Z');
mu = B(:, 1);
autoregressive_matrices{1} = B(:, 2:end);

% Sims
% (provides same result as above)
% var = rfvar3(Y',1,zeros(size(Y')),0,5,2)
% mu = var.B(:, 1);
% autoregressive_matrices{1} = var.B(:, 2:end);

%% save values
save('my_var_est.mat', 'mu', 'autoregressive_matrices');

%% Call Dynare with VAR
clearvars
clearvars -global
dynare example1_var.mod
save('nkm_var_saved_data.mat')

%% compare values
clearvars
clearvars -global
zerotol = 1e-12;

nv = load('nkm_saved_data.mat');
wv = load('nkm_var_saved_data.mat');

ridx = 3;
cidx = 2;
exo_names = nv.M_.exo_names;
endo_names = nv.M_.endo_names;

for i = 1:length(exo_names)
    figure('Name', ['Shock to ' exo_names(i)]);
    for j = 1:length(endo_names)
        subplot(ridx, cidx, j);
        hold on
        title(endo_names(j));
        plot(nv.oo_.irfs.([endo_names(j) '_' exo_names(i)]));
        plot(wv.oo_.irfs.([endo_names(j) '_' exo_names(i)]), '--');
        hold off
    end
end