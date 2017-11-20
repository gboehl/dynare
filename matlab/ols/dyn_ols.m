function ds = dyn_ols(ds, fitted_names_dict, eqtags)
% function ds = dyn_ols(ds, fitted_names_dict, eqtags)
% Run OLS on chosen model equations; unlike olseqs, allow for time t
% endogenous variables on LHS
%
% INPUTS
%   ds                [dseries]    data
%   fitted_names_dict [cell]       Nx2 cell array to be used in naming fitted
%                                  values; first column is the var name,
%                                  second column is the name of the
%                                  associated fitted value.
%   eqtags            [cellstr]    names of equation tags to estimate. If empty,
%                                  estimate all equations
%
% OUTPUTS
%   ds                [dseries]    data updated with fitted values
%
% SPECIAL REQUIREMENTS
%   none

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

assert(isdseries(ds), 'dyn_ols: the first argument must be a dseries');

jsonfile = [M_.fname '_original.json'];
if exist(jsonfile, 'file') ~= 2
    error('Could not find %s! Please use the json option (See the Dynare invocation section in the reference manual).', jsonfile);
end

%% Get Equation(s)
jsonmodel = loadjson(jsonfile);
jsonmodel = jsonmodel.model;

if nargin == 1
    [lhs, rhs, lineno, sample] = getEquationsByTags(jsonmodel);
    fitted_names_dict = {};
else
    assert(isempty(fitted_names_dict) || ...
        (iscell(fitted_names_dict) && columns(fitted_names_dict) == 2), ...
        'dyn_ols: the second argument must be an Nx2 cell array');
    if nargin == 2
        [lhs, rhs, lineno, sample] = getEquationsByTags(jsonmodel);
    else
        [lhs, rhs, lineno, sample] = getEquationsByTags(jsonmodel, 'name', eqtags);
    end
    if isempty(lhs)
        disp('dyn_ols: Nothing to estimate')
        return
    end
end

%% Estimation
M_endo_trim = cellstr(M_.endo_names);
M_exo_trim = cellstr(M_.exo_names);
M_endo_exo_names_trim = [M_endo_trim; M_exo_trim];
regex = strjoin(M_endo_exo_names_trim(:,1), '|');
mathops = '[\+\*\^\-\/\(\)]';
M_param_names_trim = cellfun(@strtrim, num2cell(M_.param_names,2), 'UniformOutput', false);
for i = 1:length(lhs)
    %% Construct regression matrices
    rhs_ = strsplit(rhs{i}, {'+','-','*','/','^','log(','exp(','(',')'});
    rhs_(cellfun(@(x) all(isstrprop(x, 'digit')), rhs_)) = [];
    vnames = setdiff(rhs_, cellstr(M_.param_names));
    if ~isempty(regexp(rhs{i}, ...
            ['(' strjoin(vnames, '\\(\\d+\\)|') '\\(\\d+\\))'], ...
            'once'))
        error(['dyn_ols: you cannot have leads in equation on line ' ...
            lineno{i} ': ' lhs{i} ' = ' rhs{i}]);
    end

    pnames = intersect(rhs_, cellstr(M_.param_names));
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
        if rhs{i}(startidx) == '*'
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

    lhssub = dseries();
    rhs_ = strsplit(rhs{i}, [splitstrings; pnames]);
    for j = 1:length(rhs_)
        rhsj = rhs_{j};
        while ~isempty(rhsj)
            minusstr = '';
            if strcmp(rhsj(1), '-') || strcmp(rhsj(1), '+')
                if length(rhsj) == 1
                    break
                end
                if strcmp(rhsj(1), '-')
                    minusstr = '-';
                end
                rhsj = rhsj(2:end);
            end
            str = getStrMoveRight(rhsj);
            if ~isempty(str)
                try
                    lhssub = [lhssub eval(regexprep([minusstr str], regex, 'ds.$&'))];
                    lhssub.rename_(lhssub{lhssub.vobs}.name{:}, [minusstr str]);
                catch
                    if ~any(strcmp(M_exo_trim, str))
                        error(['dyn_ols: problem evaluating ' minusstr str]);
                    end
                end
                rhsj = rhsj(length(str)+1:end);
            end
        end
    end
 
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
    if nargin == 3
        if iscell(eqtags)
            tagv = eqtags{i};
        else
            tagv = eqtags;
        end
    else
        tagv = ['eq_line_no_' num2str(lineno{i})];
    end
    [nobs, nvars] = size(X);
    oo_.ols.(tagv).dof = nobs - nvars;

    % Estimated Parameters
    [q, r] = qr(X, 0);
    xpxi = (r'*r)\eye(nvars);
    oo_.ols.(tagv).beta = r\(q'*Y.data);
    for j = 1:length(pnames)
        M_.params(strcmp(M_param_names_trim, pnames{j})) = oo_.ols.(tagv).beta(j);
    end

    % Yhat
    lhsrep = regexprep(lhs{i}, '[\(\)\-+\*/]', '_');
    yhatname = [lhsrep '_FIT'];
    if ~isempty(fitted_names_dict)
        idx = strcmp(fitted_names_dict(:,1), lhsrep);
        if any(idx)
            yhatname = fitted_names_dict{idx, 2};
        end
    end
    oo_.ols.(tagv).Yhat = dseries(X*oo_.ols.(tagv).beta, fp, yhatname);

    % Residuals
    oo_.ols.(tagv).resid = Y - oo_.ols.(tagv).Yhat;

    % Correct Yhat reported back to user for given
    for j = 1:lhssub.vobs
        oo_.ols.(tagv).Yhat = oo_.ols.(tagv).Yhat + lhssub{j}(fp:lp);
    end
    ds = [ds oo_.ols.(tagv).Yhat];

    %% Calculate statistics
    % Estimate for sigma^2
    SS_res = oo_.ols.(tagv).resid.data'*oo_.ols.(tagv).resid.data;
    oo_.ols.(tagv).s2 = SS_res/oo_.ols.(tagv).dof;

    % R^2
    ym = Y.data - mean(Y);
    SS_tot = ym'*ym;
    oo_.ols.(tagv).R2 = 1 - SS_res/SS_tot;

    % Adjusted R^2
    oo_.ols.(tagv).adjR2 = oo_.ols.(tagv).R2 - (1 - oo_.ols.(tagv).R2)*nvars/(oo_.ols.(tagv).dof-1);

    % Durbin-Watson
    ediff = oo_.ols.(tagv).resid.data(2:nobs) - oo_.ols.(tagv).resid.data(1:nobs-1);
    oo_.ols.(tagv).dw = (ediff'*ediff)/SS_res;

    % Standard Error
    oo_.ols.(tagv).stderr = sqrt(oo_.ols.(tagv).s2*diag(xpxi));

    % T-Stat
    oo_.ols.(tagv).tstat = oo_.ols.(tagv).beta./oo_.ols.(tagv).stderr;

    %% Print Output
    if ~options_.noprint
        title = sprintf('OLS Estimation of equation  `%s`', tagv);
        if nargin == 3
            title = [title sprintf(' [%s = %s]', 'name', tagv)];
        end

        preamble = {sprintf('Dependent Variable: %s', lhs{i}), ...
            sprintf('No. Independent Variables: %d', nvars), ...
            sprintf('Observations: %d from %s to %s\n', nobs, fp.char, lp.char)};

        afterward = {sprintf('R^2: %f', oo_.ols.(tagv).R2), ...
            sprintf('R^2 Adjusted: %f', oo_.ols.(tagv).adjR2), ...
            sprintf('s^2: %f', oo_.ols.(tagv).s2), ...
            sprintf('Durbin-Watson: %f', oo_.ols.(tagv).dw)};

        dyn_table(title, preamble, afterward, vnames, ...
            {'Coefficients','t-statistic','Std. Error'}, 4, ...
            [oo_.ols.(tagv).beta oo_.ols.(tagv).tstat oo_.ols.(tagv).stderr]);
    end
end
end
