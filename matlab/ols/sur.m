function varargout = sur(ds, param_names, eqtags, model_name)
%function varargout = sur(ds, param_names, eqtags, model_name)
% Seemingly Unrelated Regressions
%
% INPUTS
%   ds                [dseries]    data to use in estimation
%   param_names       [cellstr]    list of parameters to estimate
%   eqtags            [cellstr]    names of equation tags to estimate. If empty,
%                                  estimate all equations
%   model_name        [string]     name of model to be displayed with
%                                  report
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   dynare must have been run with the option: json=compute

% Copyright (C) 2017-2019 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global M_ oo_ options_

%% Check input argument
if nargin < 1 || nargin > 4
    error('function takes between 1 and 4 arguments');
end

if nargin < 4
    if ~isfield(oo_, 'sur')
        model_name = 'sur_model_number_1';
    else
        model_name = ['sur_model_number_' num2str(length(fieldnames(oo_.sur)) + 1)];
    end
end

if nargin < 3
    eqtags = {};
end

if nargin < 2
    param_names = {};
else
    assert(iscellstr(param_names), 'the 2nd argument must be a cellstr');
end

%% Get Equation(s)
ast = handle_constant_eqs(get_ast(eqtags));
neqs = length(ast);

%% Find parameters and variable names in equations and setup estimation matrices
[Y, lhssub, X] = common_parsing(ds, ast, true);
clear ast
nobs = Y{1}.nobs;
[Y, lhssub, X, constrained] = put_in_sur_form(Y, lhssub, X);

if nargin == 1 && size(X, 2) ~= M_.param_nbr
    warning(['Not all parameters were used in model: ' strjoin(setdiff(M_.param_names, X.name), ', ')]);
end

if ~isempty(param_names)
    newlhssub = dseries();
    nparams = length(param_names);
    pidxs = zeros(nparams, 1);
    names = X.name;
    for i = 1:nparams
        pidxs(i) = find(strcmp(param_names{i}, names));
        if isempty(pidxs(i))
            if ~isempty(eqtags)
                error(['Could not find ' param_names{i} ...
                    ' in the equations specified by ' strjoin(eqtags, ',')]);
            end
            error('Couldn''t find parameter in equations');
        end
    end
    subcols = setdiff(1:X.vobs, pidxs);
    for i = 1:length(subcols)
        newlhssub = newlhssub + M_.params(strcmp(names{subcols(i)}, M_.param_names)) * X.(names{subcols(i)});
        X = X.remove(names{subcols(i)});
    end
    Y = Y - newlhssub;
    lhssub = lhssub + newlhssub;
end

% opidxs: indexes in M_.params associated with columns of X
opidxs = zeros(X.vobs, 1);
for i = 1:X.vobs
    opidxs(i, 1) = find(strcmp(X.name{i}, M_.param_names));
end

%% Return to surgibbs if called from there
st = dbstack(1);
if strcmp(st(1).name, 'surgibbs')
    varargout{1} = nobs;
    varargout{2} = opidxs;
    varargout{3} = X{param_names{:}}.data;
    varargout{4} = Y.data;
    varargout{5} = neqs;
    return
end

% constrained_param_idxs: indexes in X.name of parameters that were constrained
constrained_param_idxs = [];
for i = 1:length(constrained)
    idx = find(strcmp(X.name, constrained{i}));
    if ~isempty(idx)
        constrained_param_idxs(end+1, 1) = idx;
    end
end

%% Estimation
oo_.sur.(model_name).dof = nobs;

% Estimated Parameters
[q, r] = qr(X.data, 0);
xpxi = (r'*r)\eye(size(X.data, 2));
resid = Y.data - X.data * (r\(q'*Y.data));
resid = reshape(resid, oo_.sur.(model_name).dof, neqs);

M_.Sigma_e = resid'*resid/oo_.sur.(model_name).dof;
kLeye = kron(chol(inv(M_.Sigma_e)), eye(oo_.sur.(model_name).dof));
[q, r] = qr(kLeye*X.data, 0);
oo_.sur.(model_name).beta = r\(q'*kLeye*Y.data);

M_.params(opidxs) = oo_.sur.(model_name).beta;

% Yhat
oo_.sur.(model_name).Yhat = X.data * oo_.sur.(model_name).beta;

% Residuals
oo_.sur.(model_name).resid = Y.data - oo_.sur.(model_name).Yhat;

% Correct Yhat reported back to user
oo_.sur.(model_name).Yhat = oo_.sur.(model_name).Yhat + lhssub;

%% Calculate statistics
% Estimate for sigma^2
SS_res = oo_.sur.(model_name).resid'*oo_.sur.(model_name).resid;
oo_.sur.(model_name).s2 = SS_res/oo_.sur.(model_name).dof;

% R^2
ym = Y.data - mean(Y.data);
SS_tot = ym'*ym;
oo_.sur.(model_name).R2 = 1 - SS_res/SS_tot;

% Adjusted R^2
oo_.sur.(model_name).adjR2 = oo_.sur.(model_name).R2 - (1 - oo_.sur.(model_name).R2)*M_.param_nbr/(oo_.sur.(model_name).dof - 1);

% Durbin-Watson
ediff = oo_.sur.(model_name).resid(2:oo_.sur.(model_name).dof) - oo_.sur.(model_name).resid(1:oo_.sur.(model_name).dof - 1);
oo_.sur.(model_name).dw = (ediff'*ediff)/SS_res;

% Standard Error
oo_.sur.(model_name).stderr = sqrt(oo_.sur.(model_name).s2*diag(xpxi));

% T-Stat
oo_.sur.(model_name).tstat = oo_.sur.(model_name).beta./oo_.sur.(model_name).stderr;

oo_.sur.(model_name).neqs = neqs;
oo_.sur.(model_name).pname = X.name;

%% Print Output
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
