function varargout = sur(ds, param_names, eqtags, model_name, noniterative, ds_range)

% Seemingly Unrelated Regressions
%
% INPUTS
%   ds                [dseries]    data to use in estimation
%   param_names       [cellstr]    list of parameters to estimate
%   eqtags            [cellstr]    names of equation tags to estimate. If empty,
%                                  estimate all equations
%   model_name        [string]     name of model to be displayed with
%                                  report
%   noniterative      [bool]       if true use noniterative estimation method
%   ds_range          [dates]      range of dates to use in estimation
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   dynare must have been run with the option: json=compute

% Copyright © 2017-2023 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

global M_ oo_ options_

%
% Check inputs
%

if nargin < 1 || nargin > 6
    error('function takes between 1 and 6 arguments');
end

if isempty(ds) || ~isdseries(ds)
    error('The first argument must be a dseries');
end

if nargin < 6
    ds_range = ds.dates;
else
    if isempty(ds_range)
        ds_range = ds.dates;
    else
        if ds_range(1) < ds.firstdate || ds_range(end) > lastdate(ds)
            error('There is a problem with the 6th argument: the date range does not correspond to that of the dseries')
        end
    end
end

if nargin < 5
    noniterative = false;
else
    if ~islogical(noniterative)
        error('the fifth argument, if passed, must be a boolean')
    end
end

if nargin < 4
    if ~isfield(oo_, 'sur')
        model_name = 'sur_model_number_1';
    else
        model_name = ['sur_model_number_' num2str(length(fieldnames(oo_.sur)) + 1)];
    end
else
    if ~isvarname(model_name)
        error('the fourth argument, if passed, must be a string')
    end
end

if nargin < 3
    eqtags = {};
end

if nargin < 2
    param_names = {};
else
    if ~isempty(param_names) && ~iscellstr(param_names)
        error('the 2nd argument must be a cellstr')
    end
end

maxit = 100;
tol = 1e-6;

%
% Get Equation(s)
%

ast = handle_constant_eqs(get_ast(eqtags));
neqs = length(ast);

%
% Find parameters and variable names in equations and setup estimation matrices
%

[Y, lhssub, X, fp, lp, residnames] = common_parsing(ds(ds_range), ast, true, param_names);
clear ast
nobs = Y{1}.nobs;
[Y, lhssub, X, constrained] = put_in_sur_form(Y, lhssub, X);

if nargin == 1 && size(X, 2) ~= M_.param_nbr
    warning(['Not all parameters were used in model: ' strjoin(setdiff(M_.param_names, X.name)', ', ')]);
end

%
% Return to surgibbs if called from there
%

st = dbstack(1);
if ~isempty(st) && strcmp(st(1).name, 'surgibbs')
    varargout{1} = nobs;
    varargout{2} = X{param_names{:}}.data;
    varargout{3} = Y.data;
    varargout{4} = neqs;
    varargout{5} = lhssub.data;
    varargout{6} = fp{1};
    varargout{7} = lp{1};
    return
end

% constrained_param_idxs: indexes in X.name of parameters that were constrained
constrained_param_idxs = NaN(length(constrained), 1);
j = 0;
for i = 1:length(constrained)
    idx = find(strcmp(X.name, constrained{i}));
    if ~isempty(idx)
        j = j+1;
        constrained_param_idxs(j, 1) = idx;
    end
end
constrained_param_idxs = constrained_param_idxs(1:j);

%
% Estimation
%

oo_.sur.(model_name).dof = nobs;

% Estimated Parameters
[q, r] = qr(X.data, 0);
xpxi = (r'*r)\eye(size(X.data, 2));
beta0 = (r\(q'*Y.data));
for i = 1:maxit
    resid = Y.data - X.data * beta0;
    resid = reshape(resid, oo_.sur.(model_name).dof, neqs);
    vcv = resid'*resid/oo_.sur.(model_name).dof;
    kLeye = kron(inv(chol(vcv))', eye(oo_.sur.(model_name).dof));
    [q, r] = qr(kLeye*X.data, 0);
    oo_.sur.(model_name).beta = r\(q'*kLeye*Y.data);
    if noniterative || max(abs(beta0 - oo_.sur.(model_name).beta)) < tol
        break
    end
    beta0 = oo_.sur.(model_name).beta;
    if i == maxit
        warning('maximum number of iterations reached')
    end
end

% Set appropriate entries in M_.Sigma_e
idxs = zeros(length(residnames), 1);
for i = 1:length(residnames)
    idxs(i) = find(strcmp(residnames{i}, M_.exo_names));
end
M_.Sigma_e(idxs, idxs) = vcv;

% opidxs: indexes in M_.params associated with columns of X
opidxs = zeros(X.vobs, 1);
for i = 1:X.vobs
    opidxs(i, 1) = find(strcmp(X.name{i}, M_.param_names));
end

% Set params
M_.params(opidxs) = oo_.sur.(model_name).beta;

% Write .inc file
write_param_init_inc_file('sur', model_name, opidxs, oo_.sur.(model_name).beta);

% Yhat
oo_.sur.(model_name).Yhat = X.data * oo_.sur.(model_name).beta;
oo_.sur.(model_name).YhatOrig = oo_.sur.(model_name).Yhat;
oo_.sur.(model_name).Yobs = Y;

% Residuals
oo_.sur.(model_name).resid = Y.data - oo_.sur.(model_name).Yhat;

% Correct Yhat reported back to user
oo_.sur.(model_name).Yhat = oo_.sur.(model_name).Yhat + lhssub;
yhatname = [model_name '_FIT'];
ds.(yhatname) = dseries(oo_.sur.(model_name).Yhat.data,  fp{1}, yhatname);
varargout{1} = ds;

%
% Calculate various statistics
%

% Estimate for sigma^2
SS_res = oo_.sur.(model_name).resid'*oo_.sur.(model_name).resid;
oo_.sur.(model_name).s2 = SS_res/oo_.sur.(model_name).dof;

% System R^2 value of McElroy (1977) - formula from Judge et al. (1986, p. 477)
oo_.sur.(model_name).R2 = 1 - (oo_.sur.(model_name).resid' * kron(inv(M_.Sigma_e(idxs,idxs)), eye(nobs)) * oo_.sur.(model_name).resid) ...
                            / (oo_.sur.(model_name).Yobs.data' * kron(inv(M_.Sigma_e(idxs,idxs)), eye(nobs)-ones(nobs,nobs)/nobs) * oo_.sur.(model_name).Yobs.data);

% Adjusted R^2
oo_.sur.(model_name).adjR2 = 1 - (1 - oo_.sur.(model_name).R2) * ((neqs*nobs-neqs)/(neqs*nobs-size(oo_.sur.(model_name).beta,1)));

% Durbin-Watson
ediff = oo_.sur.(model_name).resid(2:oo_.sur.(model_name).dof) - oo_.sur.(model_name).resid(1:oo_.sur.(model_name).dof - 1);
oo_.sur.(model_name).dw = (ediff'*ediff)/SS_res;

% Standard Error
oo_.sur.(model_name).stderr = sqrt(oo_.sur.(model_name).s2*diag(xpxi));

% T-Stat
oo_.sur.(model_name).tstat = oo_.sur.(model_name).beta./oo_.sur.(model_name).stderr;

oo_.sur.(model_name).neqs = neqs;
oo_.sur.(model_name).pname = X.name;

%
% Print Output
%

if ~options_.noprint
    preamble = {['Model name: ' model_name], ...
        sprintf('No. Equations: %d', oo_.sur.(model_name).neqs), ...
        sprintf('No. Independent Variables: %d', size(X, 2)), ...
        sprintf('Observations: %d', oo_.sur.(model_name).dof)};

    afterward = {sprintf('R^2: %f', oo_.sur.(model_name).R2), ...
        sprintf('R^2 Adjusted: %f', oo_.sur.(model_name).adjR2), ...
        sprintf('s^2: %f', oo_.sur.(model_name).s2), ...
        sprintf('Durbin-Watson: %f', oo_.sur.(model_name).dw)};

    if ~isempty(constrained_param_idxs)
        afterward = [afterward, ['Constrained parameters: ' ...
            strjoin(X.name(constrained_param_idxs), ', ')]];
    end

    dyn_table('SUR Estimation', preamble, afterward, X.name, ...
        {'Estimates','t-statistic','Std. Error'}, 4, ...
        [oo_.sur.(model_name).beta oo_.sur.(model_name).tstat oo_.sur.(model_name).stderr]);
end
end
