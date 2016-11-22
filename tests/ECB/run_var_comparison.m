clearvars
clearvars -global
!rm nkm_saved_data.mat
!rm my_var_est.mat
!rm nkm_var_saved_data.mat

%% Call Dynare without VAR
dynare nkm.mod
save('nkm_saved_data.mat')

%% Estimate VAR using simulated data
disp('VAR Estimation');
Y = oo_.endo_simul(2:3, 2:end);
Z = [ones(1, size(Y,2)); oo_.endo_simul(2:3, 1:end-1)];

% OLS
B = Y*transpose(Z)/(Z*transpose(Z));

%% save values
mu = B(:, 1);
autoregressive_matrices{1} = B(:, 2:end);
save('my_var_est.mat', 'mu', 'autoregressive_matrices');

%% Call Dynare with VAR
clearvars
clearvars -global
dynare nkm_var.mod
save('nkm_var_saved_data.mat')

%% compare values
clearvars
clearvars -global
zerotol = 1e-12;

nv = load('nkm_saved_data.mat');
wv = load('nkm_var_saved_data.mat');

assert(max(max(abs(nv.y - wv.y))) < zerotol);
assert(max(max(abs(nv.pi - wv.pi))) < zerotol);
assert(max(max(abs(nv.i - wv.i))) < zerotol);

fn = fieldnames(nv.oo_.irfs);
for i=1:length(fn)
    assert(max(max(abs(nv.oo_.irfs.(fn{i}) - (wv.oo_.irfs.(fn{i}))))) < zerotol);
end