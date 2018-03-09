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
    fitted_names_dict = {};
elseif nargin == 2
    assert(isempty(fitted_names_dict) || ...
        (iscell(fitted_names_dict) && ...
        (size(fitted_names_dict, 2) == 2 || size(fitted_names_dict, 2) == 3)), ...
        'dyn_ols: the second argument must be an Nx2 or Nx3 cell array');
elseif nargin == 3
    jsonmodel = getEquationsByTags(jsonmodel, 'name', eqtags);
end

%% Estimation
M_endo_exo_names_trim = [M_.endo_names; M_.exo_names];
[junk, idxs] = sort(cellfun(@length, M_endo_exo_names_trim), 'descend');
regex = strjoin(M_endo_exo_names_trim(idxs), '|');
mathops = '[\+\*\^\-\/\(\)]';
for i = 1:length(jsonmodel)
    %% Construct regression matrices
    rhs_ = strsplit(jsonmodel{i}.rhs, {'+','-','*','/','^','log(','exp(','(',')'});
    rhs_(cellfun(@(x) all(isstrprop(x, 'digit')), rhs_)) = [];
    vnames = setdiff(rhs_, M_.param_names);
    if ~isempty(regexp(jsonmodel{i}.rhs, ...
            ['(' strjoin(vnames, '\\(\\d+\\)|') '\\(\\d+\\))'], ...
            'once'))
        error(['dyn_ols: you cannot have leads in equation on line ' ...
            jsonmodel{i}.line ': ' jsonmodel{i}.lhs ' = ' jsonmodel{i}.rhs]);
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
        [startidx, endidx] = regexp(jsonmodel{i}.rhs, pregex, 'start', 'end');
        assert(length(startidx) == 1);
        if jsonmodel{i}.rhs(startidx) == '*' && jsonmodel{i}.rhs(endidx) == '*'
            vnamesl = getStrMoveLeft(jsonmodel{i}.rhs(1:startidx-1));
            vnamesr = getStrMoveRight(jsonmodel{i}.rhs(endidx+1:end));
            vnames{j} = [vnamesl '*' vnamesr];
            splitstrings{j} = [vnamesl '*' pnames{j} '*' vnamesr];
        elseif jsonmodel{i}.rhs(startidx) == '*'
            vnames{j} = getStrMoveLeft(jsonmodel{i}.rhs(1:startidx-1));
            splitstrings{j} = [vnames{j} '*' pnames{j}];
        elseif jsonmodel{i}.rhs(endidx) == '*'
            vnames{j} = getStrMoveRight(jsonmodel{i}.rhs(endidx+1:end));
            splitstrings{j} = [pnames{j} '*' vnames{j}];
            if jsonmodel{i}.rhs(startidx) == '-'
                vnames{j} = ['-' vnames{j}];
                splitstrings{j} = ['-' splitstrings{j}];
            end
        elseif jsonmodel{i}.rhs(startidx) == '+' ...
                || jsonmodel{i}.rhs(startidx) == '-' ...
                || jsonmodel{i}.rhs(endidx) == '+' ...
                || jsonmodel{i}.rhs(endidx) == '-'
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
            if jsonmodel{i}.rhs(startidx) == '-'
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

    lhssub = getRhsToSubFromLhs(ds, jsonmodel{i}.rhs, regex, splitstrings, pnames);
    residuals = setdiff(intersect(rhs_, M_.exo_names), ds.name);
    assert(~isempty(residuals), ['No residuals in equation ' num2str(i)]);
    assert(length(residuals) == 1, ['More than one residual in equation ' num2str(i)]);

    Y = eval(regexprep(jsonmodel{i}.lhs, regex, 'ds.$&'));
    if ~isempty(lhssub)
        Y = Y - lhssub;
    end

    fp = max(Y.firstobservedperiod, X.firstobservedperiod);
    lp = min(Y.lastobservedperiod, X.lastobservedperiod);
    if isfield(jsonmodel{i}, 'tags') ...
            && isfield(jsonmodel{i}.tags, 'sample') ...
            && ~isempty(jsonmodel{i}.tags.sample)
        colon_idx = strfind(jsonmodel{i}.tags.sample, ':');
        fsd = dates(jsonmodel{i}.tags.sample(1:colon_idx-1));
        lsd = dates(jsonmodel{i}.tags.sample(colon_idx+1:end));
        if fp > fsd
            warning(['The sample over which you want to estimate contains NaNs. '...
                'Adjusting estimation range to begin on: ' fp.char])
        else
            fp = fsd;
        end
        if lp < lsd
             warning(['The sample over which you want to estimate contains NaNs. '...
                'Adjusting estimation range to end on: ' lp.char])
        else
            lp = lsd;
        end
    end

    Y = Y(fp:lp);
    X = X(fp:lp).data;
    if ~isempty(lhssub)
        lhssub = lhssub(fp:lp);
    end

    if isfield(jsonmodel{i}, 'tags') && ...
            isfield(jsonmodel{i}.tags, 'name')
        tag = jsonmodel{i}.tags.('name');
    else
        tag = ['eq_line_no_' num2str(jsonmodel{i}.line)];
    end

    %% Estimation
    % From LeSage, James P. "Applied Econometrics using MATLAB"
    [nobs, nvars] = size(X);
    oo_.ols.(tag).dof = nobs - nvars;

    % Estimated Parameters
    [q, r] = qr(X, 0);
    xpxi = (r'*r)\eye(nvars);
    oo_.ols.(tag).beta = r\(q'*Y.data);
    oo_.ols.(tag).param_idxs = zeros(length(pnames), 1);
    for j = 1:length(pnames)
        oo_.ols.(tag).param_idxs(j) = find(strcmp(M_.param_names, pnames{j}));
        M_.params(oo_.ols.(tag).param_idxs(j)) = oo_.ols.(tag).beta(j);
    end

    % Yhat
    idx = 0;
    yhatname = [tag '_FIT'];
    if ~isempty(fitted_names_dict)
        idx = strcmp(fitted_names_dict(:,1), tag);
        if any(idx)
            yhatname = fitted_names_dict{idx, 2};
        end
    end
    oo_.ols.(tag).Yhat = dseries(X*oo_.ols.(tag).beta, fp, yhatname);

    % Residuals
    oo_.ols.(tag).resid = Y - oo_.ols.(tag).Yhat;

    % Correct Yhat reported back to user
    if ~isempty(lhssub)
        Y = Y + lhssub;
        oo_.ols.(tag).Yhat = oo_.ols.(tag).Yhat + lhssub;
    end
 
    % Apply correcting function for Yhat if it was passed
    if any(idx) ...
            && length(fitted_names_dict(idx, :)) == 3 ...
            && ~isempty(fitted_names_dict{idx, 3})
        oo_.ols.(tag).Yhat = ...
            feval(fitted_names_dict{idx, 3}, oo_.ols.(tag).Yhat);
    end
    ds = [ds oo_.ols.(tag).Yhat];

    %% Calculate statistics
    % Estimate for sigma^2
    SS_res = oo_.ols.(tag).resid.data'*oo_.ols.(tag).resid.data;
    oo_.ols.(tag).s2 = SS_res/oo_.ols.(tag).dof;

    % R^2
    ym = Y.data - mean(Y);
    SS_tot = ym'*ym;
    oo_.ols.(tag).R2 = 1 - SS_res/SS_tot;

    % Adjusted R^2
    oo_.ols.(tag).adjR2 = oo_.ols.(tag).R2 - (1 - oo_.ols.(tag).R2)*(nvars-1)/(oo_.ols.(tag).dof);

    % Durbin-Watson
    ediff = oo_.ols.(tag).resid.data(2:nobs) - oo_.ols.(tag).resid.data(1:nobs-1);
    oo_.ols.(tag).dw = (ediff'*ediff)/SS_res;

    % Standard Error
    oo_.ols.(tag).stderr = sqrt(oo_.ols.(tag).s2*diag(xpxi));

    % T-Stat
    oo_.ols.(tag).tstat = oo_.ols.(tag).beta./oo_.ols.(tag).stderr;

    %% Print Output
    if ~options_.noprint
        if nargin == 3
            title = ['OLS Estimation of equation ''' tag ''' [name = ''' tag ''']'];
        else
            title = ['OLS Estimation of equation ''' tag ''''];
        end

        preamble = {sprintf('Dependent Variable: %s', jsonmodel{i}.lhs), ...
            sprintf('No. Independent Variables: %d', nvars), ...
            sprintf('Observations: %d from %s to %s\n', nobs, fp.char, lp.char)};

        afterward = {sprintf('R^2: %f', oo_.ols.(tag).R2), ...
            sprintf('R^2 Adjusted: %f', oo_.ols.(tag).adjR2), ...
            sprintf('s^2: %f', oo_.ols.(tag).s2), ...
            sprintf('Durbin-Watson: %f', oo_.ols.(tag).dw)};

        dyn_table(title, preamble, afterward, vnames, ...
            {'Coefficients','t-statistic','Std. Error'}, 4, ...
            [oo_.ols.(tag).beta oo_.ols.(tag).tstat oo_.ols.(tag).stderr]);
    end
end
end
