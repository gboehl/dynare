function varargout = sur(ds, param_names, eqtags)
%function varargout = sur(ds, param_names, eqtags)
% Seemingly Unrelated Regressions
%
% INPUTS
%   ds                [dseries]    data to use in estimation
%   param_names       [cellstr]    list of parameters to estimate
%   eqtags            [cellstr]    names of equation tags to estimate. If empty,
%                                  estimate all equations
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   dynare must be run with the option: json=compute

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
assert(nargin >= 1 && nargin <= 3, 'sur() takes between 1 and 3 arguments');

if nargin < 3
    eqtags = {};
end

if nargin < 2
    param_names = {};
else
    assert(iscellstr(param_names), 'sur: the 2nd argument must be a cellstr');
end

%% Get Equation(s)
[ast, jsonmodel] = get_ast_jsonmodel(eqtags);
neqs = length(jsonmodel);

%% Find parameters and variable names in equations and setup estimation matrices
[Y, ~, X] = common_parsing(ds, ast, jsonmodel, true);
clear ast jsonmodel;
nobs = Y{1}.nobs;





[Y, X, constrained] = put_in_sur_form(Y, X);

if nargin == 1 && size(X, 2) ~= M_.param_nbr
    warning(['Not all parameters were used in model: ' strjoin(setdiff(M_.param_names, X.name), ', ')]);
end

if ~isempty(param_names)
    newX = dseries();
    nparams = length(param_names);
    pidxs = zeros(nparams, 1);
    for i = 1:nparams
        idx = find(strcmp(param_names{i}, X.name));
        if isempty(idx)
            if ~isempty(eqtags)
                error(['Could not find ' param_names{i} ...
                    ' in the provided equations specified by ' strjoin(eqtags, ',')]);
            end
            error('Unspecified error. Please report');
        end
        pidxs(i) = idx;
        newX = [newX X.(X.name{idx})];
    end
    subcols = setdiff(1:length(X.name), pidxs);
    for i = length(subcols):-1:1
        Y = Y - M_.params(strcmp(X.name{subcols(i)}, M_.param_names))*X.(X.name{subcols(i)});
    end
    X = newX;
end

% opidxs: indexes in M_.params associated with columns of X
opidxs = zeros(length(X.name), 1);
for i = 1:length(X.name)
    opidxs(i, 1) = find(strcmp(X.name{i}, M_.param_names));
end

%% Return to surgibbs if called from there
st = dbstack(1);
if strcmp(st(1).name, 'surgibbs')
    varargout{1} = nobs;
    varargout{2} = opidxs;
    varargout{3} = X.data;
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
oo_.sur.dof = nobs;

% Estimated Parameters
[q, r] = qr(X.data, 0);
xpxi = (r'*r)\eye(size(X.data, 2));
resid = Y.data - X.data * (r\(q'*Y.data));
resid = reshape(resid, oo_.sur.dof, neqs);

M_.Sigma_e = resid'*resid/oo_.sur.dof;
kLeye = kron(chol(inv(M_.Sigma_e)), eye(oo_.sur.dof));
[q, r] = qr(kLeye*X.data, 0);
oo_.sur.beta = r\(q'*kLeye*Y.data);

M_.params(opidxs) = oo_.sur.beta;

% Yhat
oo_.sur.Yhat = X.data * oo_.sur.beta;

% Residuals
oo_.sur.resid = Y.data - oo_.sur.Yhat;

%% Calculate statistics
% Estimate for sigma^2
SS_res = oo_.sur.resid'*oo_.sur.resid;
oo_.sur.s2 = SS_res/oo_.sur.dof;

% R^2
ym = Y.data - mean(Y.data);
SS_tot = ym'*ym;
oo_.sur.R2 = 1 - SS_res/SS_tot;

% Adjusted R^2
oo_.sur.adjR2 = oo_.sur.R2 - (1 - oo_.sur.R2)*M_.param_nbr/(oo_.sur.dof - 1);

% Durbin-Watson
ediff = oo_.sur.resid(2:oo_.sur.dof) - oo_.sur.resid(1:oo_.sur.dof - 1);
oo_.sur.dw = (ediff'*ediff)/SS_res;

% Standard Error
oo_.sur.stderr = sqrt(oo_.sur.s2*diag(xpxi));

% T-Stat
oo_.sur.tstat = oo_.sur.beta./oo_.sur.stderr;

%% Print Output
if ~options_.noprint
    preamble = {sprintf('No. Equations: %d', neqs), ...
        sprintf('No. Independent Variables: %d', size(X, 2)), ...
        sprintf('Observations: %d', oo_.sur.dof)};

    afterward = {sprintf('R^2: %f', oo_.sur.R2), ...
        sprintf('R^2 Adjusted: %f', oo_.sur.adjR2), ...
        sprintf('s^2: %f', oo_.sur.s2), ...
        sprintf('Durbin-Watson: %f', oo_.sur.dw)};

    if ~isempty(constrained_param_idxs)
        afterward = [afterward, ['Constrained parameters: ' ...
            strjoin(X.name(constrained_param_idxs), ', ')]];
    end

    dyn_table('SUR Estimation', preamble, afterward, X.name, ...
        {'Estimates','t-statistic','Std. Error'}, 4, ...
        [oo_.sur.beta oo_.sur.tstat oo_.sur.stderr]);
end
end
