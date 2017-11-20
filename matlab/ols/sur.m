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
%   dynare must be run with the option: json=parse

% Copyright (C) 2017 Dynare Team
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
    error('Could not find %s! Please use the json=parse option (See the Dynare invocation section in the reference manual).', jsonfile);
end

jsonmodel = loadjson(jsonfile);
jsonmodel = jsonmodel.model;
[lhs, rhs, lineno] = getEquationsByTags(jsonmodel);

%% Find parameters and variable names in equations and setup estimation matrices
M_exo_names_trim = cellstr(M_.exo_names);
M_endo_exo_names_trim = [cellstr(M_.endo_names); M_exo_names_trim];
M_param_names_trim = cellstr(M_.param_names);
regex = strjoin(M_endo_exo_names_trim(:,1), '|');
mathops = '[\+\*\^\-\/]';
params = cell(length(rhs),1);
vars = cell(length(rhs),1);
Y = [];
X = [];
startidxs = zeros(length(lhs), 1);
startdates = cell(length(lhs), 1);
enddates = cell(length(lhs), 1);
residnames = cell(length(lhs), 1);
pidxs = zeros(M_.param_nbr, 1);
pidx = 0;
vnamesall = {};
for i = 1:length(lhs)
    rhs_ = strsplit(rhs{i}, {'+','-','*','/','^','log(','ln(','log10(','exp(','(',')','diff('});
    rhs_(cellfun(@(x) all(isstrprop(x, 'digit')), rhs_)) = [];
    vnames = setdiff(rhs_, M_param_names_trim);
    if ~isempty(regexp(rhs{i}, ...
            ['(' strjoin(vnames, '\\(\\d+\\)|') '\\(\\d+\\))'], ...
            'once'))
        error(['sur1: you cannot have leads in equation on line ' ...
            lineno{i} ': ' lhs{i} ' = ' rhs{i}]);
    end

    % Find parameters and associated variables
    pnames = intersect(rhs_, M_param_names_trim);
    vnames = cell(1, length(pnames));
    xjdata = dseries;
    for j = 1:length(pnames)
        pidx = pidx + 1;
        pidxs(pidx, 1) = find(strcmp(pnames{j}, M_param_names_trim));
        createdvar = false;
        pregex = [...
            mathops pnames{j} mathops ...
            '|^' pnames{j} mathops ...
            '|' mathops pnames{j} '$' ...
            ];
        [startidx, endidx] = regexp(rhs{i}, pregex, 'start', 'end');
        assert(length(startidx) == 1);
        if rhs{i}(startidx) == '*'
            vnames{j} = getStrMoveLeft(rhs{i}(1:startidx-1));
        elseif rhs{i}(endidx) == '*'
            vnames{j} = getStrMoveRight(rhs{i}(endidx+1:end));
        elseif rhs{i}(startidx) == '+' ...
                || rhs{i}(startidx) == '-' ...
                || rhs{i}(endidx) == '+' ...
                || rhs{i}(endidx) == '-'
            % intercept
            createdvar = true;
            if any(strcmp(M_endo_exo_names_trim, 'intercept'))
                [~, vnames{j}] = fileparts(tempname);
                vnames{j} = ['intercept_' vnames{j}];
                assert(~any(strcmp(M_endo_exo_names_trim, vnames{j})));
            else
                vnames{j} = 'intercept';
            end
        else
            error('sur1: Shouldn''t arrive here');
        end
        if createdvar
            xjdatatmp = dseries(ones(ds.nobs, 1), ds.firstdate, vnames{j});
        else
            xjdatatmp = eval(regexprep(vnames{j}, regex, 'ds.$&'));
            xjdatatmp.rename_(vnames{j});
        end
        xjdatatmp.rename_(num2str(j));
        xjdata = [xjdata xjdatatmp];
    end

    residuals = intersect(rhs_, cellstr(M_.exo_names));
    for j = 1:length(residuals)
        if any(strcmp(residuals{j}, vnames))
            residuals{j} = [];
        end
    end
    idx = ~cellfun(@isempty, residuals);
    assert(sum(idx) == 1, ['More than one residual in equation ' num2str(i)]);
    residnames{i} = residuals{idx};

    params{i} = pnames;
    vars{i} = vnames;

    ydata = eval(regexprep(lhs{i}, regex, 'ds.$&'));

    fp = max(ydata.firstobservedperiod, xjdata.firstobservedperiod);
    lp = min(ydata.lastobservedperiod, xjdata.lastobservedperiod);

    startidxs(i) = length(Y) + 1;
    startdates{i} = fp;
    enddates{i} = lp;
    Y(startidxs(i):startidxs(i)+lp-fp, 1) = ydata(fp:lp).data;
    X(startidxs(i):startidxs(i)+lp-fp, end+1:end+size(xjdata(fp:lp).data,2)) = xjdata(fp:lp).data;
end

assert(size(X, 2) == M_.param_nbr, 'Not all parameters were used in model');

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
Y = newY;
X = newX;

%% Estimation
% Estimated Parameters
oo_.sur.dof = length(maxfp:minlp);
[q, r] = qr(X, 0);
xpxi = (r'*r)\eye(M_.param_nbr);
resid = Y - X * (r\(q'*Y));
resid = reshape(resid, oo_.sur.dof, length(lhs));

M_.Sigma_e = resid'*resid/oo_.sur.dof;
kLeye = kron(inv(M_.Sigma_e), eye(oo_.sur.dof));
[q, r] = qr(kLeye*X, 0);
oo_.sur.beta = r\(q'*kLeye*Y);
M_.params(pidxs, 1) = oo_.sur.beta;

% Yhat
oo_.sur.Yhat = X * oo_.sur.beta;

% Residuals
oo_.sur.resid = Y - oo_.sur.Yhat;

%% Return to surgibbs if called from there
st = dbstack(1);
if strcmp(st(1).name, 'surgibbs')
    varargout{1} = oo_.sur.dof;
    varargout{2} = size(X, 2);
    varargout{3} = M_.param_names(pidxs, :);
    varargout{4} = oo_.sur.beta;
    varargout{5} = X;
    varargout{6} = Y;
    varargout{7} = length(lhs);
    return
end

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

    dyn_table('SUR Estimation', preamble, afterward, [vars{:}], ...
        {'Coefficients','t-statistic','Std. Error'}, 4, ...
        [oo_.sur.beta oo_.sur.tstat oo_.sur.stderr]);
end
end
