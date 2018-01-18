function ds = dyn_ols(ds, fitted_names_dict, eqtags)
% function ds = dyn_ols(ds, fitted_names_dict, eqtags)
% Run OLS on chosen model equations; unlike olseqs, allow for time t
% endogenous variables on LHS
%
% INPUTS
%   ds                [dseries]    data
%   fitted_names_dict [cell]       Nx2 or Nx3 cell array to be used in naming fitted
%                                  values; first column is the equation tag,
%                                  second column is the name of the
%                                  associated fitted value, third column
%                                  (if it exists) is the function name of
%                                  the transformation to perform on the
%                                  fitted value.
%   eqtags            [cellstr]    names of equation tags to estimate. If empty,
%                                  estimate all equations
%
% OUTPUTS
%   ds                [dseries]    data updated with fitted values
%
% SPECIAL REQUIREMENTS
%   none

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

assert(nargin >= 1 && nargin <= 3, 'dyn_ols: takes between 1 and 3 arguments');
assert(isdseries(ds), 'dyn_ols: the first argument must be a dseries');

jsonfile = [M_.fname '_original.json'];
if exist(jsonfile, 'file') ~= 2
    error('Could not find %s! Please use the json=compute option (See the Dynare invocation section in the reference manual).', jsonfile);
end

%% Get Equation(s)
jsonmodel = loadjson(jsonfile);
jsonmodel = jsonmodel.model;

if nargin == 1
    [lhs, rhs, lineno, sample, tags] = getEquationsByTags(jsonmodel);
    fitted_names_dict = {};
else
    assert(isempty(fitted_names_dict) || ...
        (iscell(fitted_names_dict) && ...
        (size(fitted_names_dict, 2) == 2 || size(fitted_names_dict, 2) == 3)), ...
        'dyn_ols: the second argument must be an Nx2 or Nx3 cell array');
    if nargin == 2
        [lhs, rhs, lineno, sample, tags] = getEquationsByTags(jsonmodel);
    else
        [lhs, rhs, lineno, sample, tags] = getEquationsByTags(jsonmodel, 'name', eqtags);
    end
    if isempty(lhs)
        disp('dyn_ols: Nothing to estimate')
        return
    end
end

%% Estimation
M_endo_exo_names_trim = [M_.endo_names; M_.exo_names];
regex = strjoin(M_endo_exo_names_trim(:,1), '|');
mathops = '[\+\*\^\-\/\(\)]';
for i = 1:length(lhs)
    %% Construct regression matrices
    rhs_ = strsplit(rhs{i}, {'+','-','*','/','^','log(','exp(','(',')'});
    rhs_(cellfun(@(x) all(isstrprop(x, 'digit')), rhs_)) = [];
    vnames = setdiff(rhs_, M_.param_names);
    if ~isempty(regexp(rhs{i}, ...
            ['(' strjoin(vnames, '\\(\\d+\\)|') '\\(\\d+\\))'], ...
            'once'))
        error(['dyn_ols: you cannot have leads in equation on line ' ...
            lineno{i} ': ' lhs{i} ' = ' rhs{i}]);
    end

    pnames = intersect(rhs_, M_.param_names);
    vnames = cell(1, length(pnames));
    splitstrings = cell(length(pnames), 1);
    X = dseries();
    for j = 1:length(pnames)
        createdvar = false;
        pregex = [...
            mathops pnames{j} mathops ...
            '|^' pnames{j} mathops ...
            '|' mathops pnames{j} '$' ...
            ];
        [startidx, endidx] = regexp(rhs{i}, pregex, 'start', 'end');
        assert(length(startidx) == 1);
        if rhs{i}(startidx) == '*' && rhs{i}(endidx) == '*'
            vnamesl = getStrMoveLeft(rhs{i}(1:startidx-1));
            vnamesr = getStrMoveRight(rhs{i}(endidx+1:end));
            vnames{j} = [vnamesl '*' vnamesr];
            splitstrings{j} = [vnamesl '*' pnames{j} '*' vnamesr];
        elseif rhs{i}(startidx) == '*'
            vnames{j} = getStrMoveLeft(rhs{i}(1:startidx-1));
            splitstrings{j} = [vnames{j} '*' pnames{j}];
        elseif rhs{i}(endidx) == '*'
            vnames{j} = getStrMoveRight(rhs{i}(endidx+1:end));
            splitstrings{j} = [pnames{j} '*' vnames{j}];
            if rhs{i}(startidx) == '-'
                vnames{j} = ['-' vnames{j}];
                splitstrings{j} = ['-' splitstrings{j}];
            end
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
            splitstrings{j} = vnames{j};
        else
            error('dyn_ols: Shouldn''t arrive here');
        end
        if createdvar
            if rhs{i}(startidx) == '-'
                Xtmp = dseries(-ones(ds.nobs, 1), ds.firstdate, vnames{j});
            else
                Xtmp = dseries(ones(ds.nobs, 1), ds.firstdate, vnames{j});
            end
        else
            Xtmp = eval(regexprep(vnames{j}, regex, 'ds.$&'));
            Xtmp.rename_(vnames{j});
        end
        X = [X Xtmp];
    end

    lhssub = getRhsToSubFromLhs(ds, rhs{i}, regex, [splitstrings; pnames]);
    residuals = setdiff(intersect(rhs_, M_.exo_names), ds.name);
    assert(~isempty(residuals), ['No residuals in equation ' num2str(i)]);
    assert(length(residuals) == 1, ['More than one residual in equation ' num2str(i)]);

    Y = eval(regexprep(lhs{i}, regex, 'ds.$&'));
    for j = 1:lhssub.vobs
        Y = Y - lhssub{j};
    end

    fp = max(Y.firstobservedperiod, X.firstobservedperiod);
    lp = min(Y.lastobservedperiod, X.lastobservedperiod);
    if ~isempty(sample{i})
        if fp > sample{i}(1) || lp < sample{i}(end)
            warning(['The sample over which you want to estimate contains NaNs. '...
                'Adjusting estimation range to be: ' fp.char ' to ' lp.char])
        else
            fp = sample{i}(1);
            lp = sample{i}(end);
        end
    end

    Y = Y(fp:lp);
    X = X(fp:lp).data;

    %% Estimation
    % From LeSage, James P. "Applied Econometrics using MATLAB"
    [nobs, nvars] = size(X);
    oo_.ols.(tags{i}).dof = nobs - nvars;

    % Estimated Parameters
    [q, r] = qr(X, 0);
    xpxi = (r'*r)\eye(nvars);
    oo_.ols.(tags{i}).beta = r\(q'*Y.data);
    for j = 1:length(pnames)
        M_.params(strcmp(M_.param_names, pnames{j})) = oo_.ols.(tags{i}).beta(j);
    end

    % Yhat
    idx = 0;
    yhatname = [tags{i} '_FIT'];
    if ~isempty(fitted_names_dict)
        idx = strcmp(fitted_names_dict(:,1), tags{i});
        if any(idx)
            yhatname = fitted_names_dict{idx, 2};
        end
    end
    oo_.ols.(tags{i}).Yhat = dseries(X*oo_.ols.(tags{i}).beta, fp, yhatname);
    if any(idx) ...
            && length(fitted_names_dict(idx, :)) == 3 ...
            && ~isempty(fitted_names_dict{idx, 3})
        oo_.ols.(tags{i}).Yhat = ...
            eval([fitted_names_dict{idx, 3} '(oo_.ols.(tags{' num2str(i) '}).Yhat)']);
    end

    % Residuals
    oo_.ols.(tags{i}).resid = Y - oo_.ols.(tags{i}).Yhat;

    % Correct Yhat reported back to user for given
    for j = 1:lhssub.vobs
        oo_.ols.(tags{i}).Yhat = oo_.ols.(tags{i}).Yhat + lhssub{j}(fp:lp);
    end
    ds = [ds oo_.ols.(tags{i}).Yhat];

    %% Calculate statistics
    % Estimate for sigma^2
    SS_res = oo_.ols.(tags{i}).resid.data'*oo_.ols.(tags{i}).resid.data;
    oo_.ols.(tags{i}).s2 = SS_res/oo_.ols.(tags{i}).dof;

    % R^2
    ym = Y.data - mean(Y);
    SS_tot = ym'*ym;
    oo_.ols.(tags{i}).R2 = 1 - SS_res/SS_tot;

    % Adjusted R^2
    oo_.ols.(tags{i}).adjR2 = oo_.ols.(tags{i}).R2 - (1 - oo_.ols.(tags{i}).R2)*nvars/(oo_.ols.(tags{i}).dof-1);

    % Durbin-Watson
    ediff = oo_.ols.(tags{i}).resid.data(2:nobs) - oo_.ols.(tags{i}).resid.data(1:nobs-1);
    oo_.ols.(tags{i}).dw = (ediff'*ediff)/SS_res;

    % Standard Error
    oo_.ols.(tags{i}).stderr = sqrt(oo_.ols.(tags{i}).s2*diag(xpxi));

    % T-Stat
    oo_.ols.(tags{i}).tstat = oo_.ols.(tags{i}).beta./oo_.ols.(tags{i}).stderr;

    %% Print Output
    if ~options_.noprint
        if nargin == 3
            title = ['OLS Estimation of equation ''' tags{i} ''' [name = ''' tags{i} ''']'];
        else
            title = ['OLS Estimation of equation ''' tags{i} ''''];
        end

        preamble = {sprintf('Dependent Variable: %s', lhs{i}), ...
            sprintf('No. Independent Variables: %d', nvars), ...
            sprintf('Observations: %d from %s to %s\n', nobs, fp.char, lp.char)};

        afterward = {sprintf('R^2: %f', oo_.ols.(tags{i}).R2), ...
            sprintf('R^2 Adjusted: %f', oo_.ols.(tags{i}).adjR2), ...
            sprintf('s^2: %f', oo_.ols.(tags{i}).s2), ...
            sprintf('Durbin-Watson: %f', oo_.ols.(tags{i}).dw)};

        dyn_table(title, preamble, afterward, vnames, ...
            {'Coefficients','t-statistic','Std. Error'}, 4, ...
            [oo_.ols.(tags{i}).beta oo_.ols.(tags{i}).tstat oo_.ols.(tags{i}).stderr]);
    end
end
end
