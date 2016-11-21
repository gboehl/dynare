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
fprintf(fid, '%%function ret = %s(y, h)\n', basename);
fprintf(fid, '%% Calculates the %d-step-ahead forecast from the VAR model %s\n', max(horizon), var_model_name);
fprintf(fid, '%%\n%% Created automatically by Dynare\n%%\n\n');
fprintf(fid, '%%%% Construct ret\n');
fprintf(fid, 'assert(length(y) == %d);\n', sum(sum(M_.lead_lag_incidence ~= 0)));

endo_names = cellstr(M_.endo_names);
yidx = zeros(size(endo_names));
for i=1:length(M_.var.(var_model_name).var_list_)
    yidx = yidx | strcmp(strtrim(M_.var.(var_model_name).var_list_(i,:)), endo_names);
end
fprintf(fid, 'y = y(logical([');
fprintf(fid, '%d;', yidx);
fprintf(fid, ']), :);\n');

lm = length(mu);
lc = length(autoregressive_matrices);
assert(lc == M_.var.(var_model_name).order);
fprintf(fid, 'if size(y, 1) ~= %d || size(y, 2) ~= %d\n', lm, M_.var.(var_model_name).order);
fprintf(fid, '    error(''The dimensions of y are not correct. It should be an nvars x order matrix'');\n');
fprintf(fid, 'end\n');

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
if M_.var.(var_model_name).order > 1
    mu = [mu; zeros(lm*M_.var.(var_model_name).order-lm, 1)];
end
fprintf(fid, '%%%% Calculate h-step-ahead forecast\n');
fprintf(fid, 'retidx = 1;\n');
fprintf(fid, 'ret = zeros(%d, %d);\n', length(horizon), lm);
fprintf(fid, 'for i=1:%d\n', max(horizon));
fprintf(fid, '    y = [');
fprintf(fid, [repmat(' %f ', 1, size(mu, 2)) ';'], mu');
fprintf(fid, '] + [');
fprintf(fid, [repmat(' %f ', 1, size(A, 2)) ';'], A');
fprintf(fid, ']*y(:);\n');
if length(horizon) == 1
    fprintf(fid, '    if %d == i\n', horizon);
else
    fprintf(fid, '    if any([');
    fprintf(fid, '%d ', horizon);
    fprintf(fid, '] == i)\n');
end
fprintf(fid, '        ret(retidx, :) = y(1:%d);\n', lm);
fprintf(fid, '        retidx = retidx + 1;\n');
fprintf(fid, '    end\n');
fprintf(fid, 'end\n');

% retidx = find(strcmp(dwrt, endo_names) & yidx == 1);
% assert(~isempty(retidx))
% if retidx == 1
%     fprintf(fid, 'y = y(1);\n');
% else
%     fprintf(fid, 'y = y(%d);\n', sum(yidx(1:retidx-1))+1);
% end

%% close file
fprintf(fid, 'end\n');
fclose(fid);
end