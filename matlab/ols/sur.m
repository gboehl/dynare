function varargout = sur(ds)
% function varargout = sur(ds)
% Seemingly Unrelated Regressions
%
% INPUTS
%   ds                  [dseries]  data to use in estimation
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   dynare must be run with the option: json=compute

% Copyright (C) 2017-2018 Dynare Team
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
assert(~isempty(ds) && isdseries(ds), 'The first argument must be a dseries');

%% Read JSON
jsonfile = [M_.fname '_original.json'];
if exist(jsonfile, 'file') ~= 2
    error('Could not find %s! Please use the json=compute option (See the Dynare invocation section in the reference manual).', jsonfile);
end

jsonmodel = loadjson(jsonfile);
jsonmodel = jsonmodel.model;
[lhs, rhs, lineno] = getEquationsByTags(jsonmodel);

%% Find parameters and variable names in equations and setup estimation matrices
M_exo_names_trim = cellstr(M_.exo_names);
M_param_names_trim = cellstr(M_.param_names);
[X, Y, startdates, enddates, startidxs, residnames, pbeta, vars, pidxs, surconstrainedparams] = ...
    pooled_sur_common(ds, lhs, rhs, lineno, M_exo_names_trim, M_param_names_trim);

if size(X, 2) ~= M_.param_nbr
    warning(['Not all parameters were used in model: ' ...
        sprintf('%s', strjoin(setdiff(M_param_names_trim, pbeta), ', '))]);
end

%% Force equations to have the same sample range
maxfp = max([startdates{:}]);
minlp = min([enddates{:}]);
nobs = minlp - maxfp;
newY = zeros(nobs*length(lhs), 1);
newX = zeros(nobs*length(lhs), columns(X));
lastidx = 1;
for i = 1:length(lhs)
    if i == length(lhs)
        yds = dseries(Y(startidxs(i):end), startdates{i});
        xds = dseries(X(startidxs(i):end, :), startdates{i});
    else
        yds = dseries(Y(startidxs(i):startidxs(i+1)-1), startdates{i});
        xds = dseries(X(startidxs(i):startidxs(i+1)-1, :), startdates{i});
    end
    newY(lastidx:lastidx + nobs, 1) = yds(maxfp:minlp).data;
    newX(lastidx:lastidx + nobs, :) = xds(maxfp:minlp, :).data;
    if i ~= length(lhs)
        lastidx = lastidx + nobs + 1;
    end
end

%% Return to surgibbs if called from there
st = dbstack(1);
if strcmp(st(1).name, 'surgibbs')
    varargout{1} = length(maxfp:minlp); %dof
    varargout{2} = pidxs;
    varargout{3} = newX;
    varargout{4} = newY;
    varargout{5} = length(lhs);
    return
end

Y = newY;
X = newX;
oo_.sur.dof = length(maxfp:minlp);

%% Estimation
% Estimated Parameters
[q, r] = qr(X, 0);
xpxi = (r'*r)\eye(size(X, 2));
resid = Y - X * (r\(q'*Y));
resid = reshape(resid, oo_.sur.dof, length(lhs));

M_.Sigma_e = resid'*resid/oo_.sur.dof;
kLeye = kron(chol(inv(M_.Sigma_e)), eye(oo_.sur.dof));
[q, r] = qr(kLeye*X, 0);
oo_.sur.beta = r\(q'*kLeye*Y);

M_.params(pidxs, 1) = oo_.sur.beta;

% Yhat
oo_.sur.Yhat = X * oo_.sur.beta;

% Residuals
oo_.sur.resid = Y - oo_.sur.Yhat;

%% Calculate statistics
% Estimate for sigma^2
SS_res = oo_.sur.resid'*oo_.sur.resid;
oo_.sur.s2 = SS_res/oo_.sur.dof;

% R^2
ym = Y - mean(Y);
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
    preamble = {sprintf('Dependent Variable: %s', lhs{i}), ...
        sprintf('No. Independent Variables: %d', M_.param_nbr), ...
        sprintf('Observations: %d', oo_.sur.dof)};

    afterward = {sprintf('R^2: %f', oo_.sur.R2), ...
        sprintf('R^2 Adjusted: %f', oo_.sur.adjR2), ...
        sprintf('s^2: %f', oo_.sur.s2), ...
        sprintf('Durbin-Watson: %f', oo_.sur.dw)};

    if ~isempty(surconstrainedparams)
        afterward = [afterward, ...
            sprintf('Constrained parameters: %s', ...
            strjoin(pbeta(surconstrainedparams), ', '))];
    end

    dyn_table('SUR Estimation', preamble, afterward, [vars{:}], ...
        {'Coefficients','t-statistic','Std. Error'}, 4, ...
        [oo_.sur.beta oo_.sur.tstat oo_.sur.stderr]);
end
end
