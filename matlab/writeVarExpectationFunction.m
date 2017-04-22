function writeVarExpectationFunction(var_model_name, horizon)
%function writeVarExpectationFunction(model_name)

%%
global M_;

%% open file
basename = ['var_forecast_' var_model_name];
fid = fopen([basename '.m'], 'w');
if fid == -1
    error(['Could not open ' basename '.m for writing']);
end

%% load .mat file
load(var_model_name, 'autoregressive_matrices', 'mu');
if ~exist('autoregressive_matrices', 'var') || ~exist('mu', 'var')
    error([var_model_name '.mat : must contain the variables autoregressive_matrices and mu']);
end

%%
fprintf(fid, 'function ret = %s(y)\n', basename);
fprintf(fid, '%%function ret = %s(y)\n', basename);
fprintf(fid, '%% Calculates the %d-step-ahead forecast from the VAR model %s\n', max(horizon), var_model_name);
fprintf(fid, '%%\n%% Created automatically by Dynare on %s\n%%\n\n', datestr(now));
fprintf(fid, '%%%% Construct y\n');
fprintf(fid, 'assert(length(y) == %d);\n', sum(sum(M_.lead_lag_incidence ~= 0)));

endo_names = cellstr(M_.endo_names);
nvars = size(M_.var.(var_model_name).var_list_,1);
var_model_order = M_.var.(var_model_name).order;
yidx = zeros(nvars, min(var_model_order, 2));
% first for order <= 2, drawing variables directly from their endo_names
for i=1:min(var_model_order, 2)
    if mod(i, 2) == 0
        ridx = 1;
    else
        ridx = 2;
    end
    for j=1:nvars
        yidx(j, i) = M_.lead_lag_incidence(ridx, strcmp(strtrim(M_.var.(var_model_name).var_list_(j,:)), endo_names)');
    end
end
yidx = yidx(:);

% then for order > 2
if var_model_order > 2
    y1idx = zeros((var_model_order - 2)*nvars, var_model_order - 2);
    for i=3:var_model_order
        for j=1:nvars
            varidx = [M_.aux_vars.orig_index] == find(strcmp(strtrim(M_.var.(var_model_name).var_list_(j,:)), endo_names)) ...
                & [M_.aux_vars.orig_lead_lag] == -i;
            cidx = [M_.aux_vars.endo_index];
            cidx = cidx(varidx);
            y1idx(j, i-2) = M_.lead_lag_incidence(2, cidx);
        end
    end
    yidx = [yidx ; y1idx(:)];
end
fprintf(fid, 'y = y([');
fprintf(fid, '%d ', yidx);
fprintf(fid, ']);\n');

lm = length(mu);
lc = length(autoregressive_matrices);
assert(lc == var_model_order);

A = zeros(lm*lc, lm*lc);
for i=1:lc
    if any([lm lm] ~= size(autoregressive_matrices{i}))
        error(['The dimensions of mu and autoregressive_matrices for ' var_model_name ' are off']);
    end
    col = lm*(i-1)+1:lm*i;
    A(1:lm, col) = autoregressive_matrices{i};
    if i ~= lc
        A(lm*i+1:lm*i+lm, col) = eye(lm, lm);
    end
end
if var_model_order > 1
    mu = [mu; zeros(lm*var_model_order-lm, 1)];
end
fprintf(fid, '\n%%%% Calculate %d-step-ahead forecast for VAR(%d) written as VAR(1)\n', max(horizon), var_model_order);
fprintf(fid, '%%  Follows Lütkepohl (2005) pg 15 & 34\n');
if max(horizon) == 1
    printInsideOfLoop(fid, mu, A, false);
    fprintf(fid, 'ret(1, :) = y(1:%d);\n', lm);
else
    fprintf(fid, 'retidx = 1;\n');
    fprintf(fid, 'ret = zeros(%d, %d);\n', length(horizon), lm);
    fprintf(fid, 'for i=1:%d\n', max(horizon));
    printInsideOfLoop(fid, mu, A, true);
    if length(horizon) == 1
        fprintf(fid, '    if %d == i\n', horizon);
    else
        fprintf(fid, '    if any([');
        fprintf(fid, '%d ', horizon);
        fprintf(fid, '] == i)\n');
    end
    fprintf(fid, '        ret(retidx, :) = y(1:%d);\n', lm);
    fprintf(fid, '        retidx = retidx + 1;\n');
%    fprintf(fid, '    ret([');
%    fprintf(fid, '%d ', horizon);
%    fprintf(fid, '] == i, :) = y(1:%d);\n', lm);
    fprintf(fid, '    end\n');
    fprintf(fid, 'end\n');
end

%% close file
fprintf(fid, 'end\n');
fclose(fid);
end

function printInsideOfLoop(fid, mu, A, inloop)
if inloop
    fs = '    ';
    ns = '       ';
    spaces = '        ';
else
    fs = '';
    ns = '   ';
    spaces = '    ';
end
fprintf(fid, '%sy = ...\n%s[ ... %% intercept\n%s', fs, spaces, ns);
fprintf(fid, [repmat('% f ', 1, size(mu, 2)) '; ...\n' ns], mu');
fprintf(fid, ' ] + ...\n%s[ ... %% autoregressive matrices\n%s', spaces, ns);
fprintf(fid, [repmat('% f ', 1, size(A, 2)) '; ...\n' ns], A');
fprintf(fid, ' ] * y;\n');
end