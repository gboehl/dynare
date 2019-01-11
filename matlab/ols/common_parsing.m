function [Y, lhssub, X, startdates, enddates, startidxs, residnames, pbeta, vars, surpidxs, surconestrainedparams] = common_parsing(ds, ast, jsonmodel, overlapping_dates)
%function [Y, lhssub, X, startdates, enddates, startidxs, residnames, pbeta, vars, surpidxs, surconstrainedparams] = common_parsing(ds, ast, jsonmodel, overlapping_dates)
%
% Code common to sur.m and pooled_ols.m
%
% INPUTS
%   ds                   [dseries]     dataset
%   ast                  [cell array]  JSON representation of abstract syntax tree
%   jsonmodel            [cell array]  JSON representation of model block
%   overlapping_dates    [boolean]     if true, dates are same across
%                                      equations
%
% OUTPUTS
%   Y                    [cell array]  dependent variables
%   lhssub               [cell array]  RHS to subtract from Y
%   X                    [cell array]  regressors
%   startdates           [cell array]  first observed period for each
%                                      equation
%   enddates             [cell array]  last observed period for each
%                                      equation
%   startidxs            [vector]      rows corresponding to each
%                                      equation's observations
%   residnames           [cell array]  name of residual in each equation
%   pbeta                [cell array]  parameter names corresponding to
%                                      columns of X
%   vars                 [cell array]  variable names corresponding to
%                                      parameters
%   surpidxs             [vector]      indexes in M_.params associated with
%                                      columns of X
%   surconstrainedparams [vector]      indexes of parameters that were
%                                      constrained
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2019 Dynare Team
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


%% Initialize variables
Y = cell(length(ast), 1);
lhssub = cell(length(ast), 1);
X = cell(length(ast), 1);
startdates = cell(length(ast), 1);
enddates = cell(length(ast), 1);

%% Loop over equations
neqs = length(ast);
for i = 1:neqs
    %% Parse equation i
    [Y{i}, lhssub{i}, X{i}] = parse_ols_style_equation(ds, ast{i}, jsonmodel{i}.line);

    %% Set start and end dates
    if overlapping_dates
        [startdates{i}, enddates{i}] = get_ols_start_end_dates(Y{i}, lhssub{i}, X{i}, jsonmodel{i});
    end
end

if overlapping_dates
    maxfp = max([startdates{:}]);
    minlp = min([enddates{:}]);
    for i = 1:neqs
        Y{i} = Y{i}(maxfp:minlp);
        X{i} = X{i}(maxfp:minlp);
        if ~isempty(lhssub{i})
            lhssub{i} = lhssub{i}(maxfp:minlp);
        end
    end
end



return
%%


M_endo_exo_names_trim = [M_.endo_names; M_.exo_names];
[junk, idxs] = sort(cellfun(@length, M_endo_exo_names_trim), 'descend');
regex = strjoin(M_endo_exo_names_trim(idxs), '|');
mathops = '[\+\*\^\-\/]';
params = cell(length(jsonmodel),1);
vars = cell(length(jsonmodel),1);
pbeta = {};
Y = [];
X = [];
startidxs = zeros(length(jsonmodel), 1);
startdates = cell(length(jsonmodel), 1);
enddates = cell(length(jsonmodel), 1);
residnames = cell(length(jsonmodel), 1);
surpidxs = zeros(M_.param_nbr, 1); 
surpidx = 0; 
surconstrainedparams = [];
for i = 1:length(jsonmodel)
    rhs_ = strsplit(jsonmodel{i}.rhs, {'+','-','*','/','^','log(','ln(','log10(','exp(','(',')','diff('});
    rhs_(cellfun(@(x) all(isstrprop(x, 'digit')), rhs_)) = [];
    vnames = setdiff(rhs_, M_.param_names);
    if ~isempty(regexp(jsonmodel{i}.rhs, ...
            ['(' strjoin(vnames, '\\(\\d+\\)|') '\\(\\d+\\))'], 'once'))
        error(['you cannot have leads in equation on line ' ...
            jsonmodel{i}.line ': ' jsonmodel{i}.lhs ' = ' jsonmodel{i}.rhs]);
    end

    % Find parameters and associated variables
    pnames = intersect(rhs_, M_.param_names);
    pidxs = zeros(length(pnames), 1);
    vnames = cell(1, length(pnames));
    splitstrings = cell(length(pnames), 1);
    xjdata = dseries;
    dropvname = true(1,length(pnames));
    for j = 1:length(pnames)
        createdvar = false;
        idx = find(strcmp(pbeta, pnames{j}));
        if isempty(idx)
            pbeta = [pbeta; pnames{j}];
            pidxs(j) = length(pbeta);
            surpidx = surpidx + 1;
            surpidxs(surpidx, 1) = find(strcmp(pnames{j}, M_.param_names));
        else
            pidxs(j) = idx;
            surconstrainedparams = [surconstrainedparams idx];
            dropvname(j) = false;
        end

        pregex = [...
            mathops pnames{j} mathops ...
            '|^' pnames{j} mathops ...
            '|' mathops pnames{j} '$' ...
            ];
        [startidx, endidx] = regexp(jsonmodel{i}.rhs, pregex, 'start', 'end');
        assert(length(startidx) == 1);
        if jsonmodel{i}.rhs(startidx) == '*' && jsonmodel{i}.rhs(endidx) == '*'
            vnames{j} = [getStrMoveLeft(jsonmodel{i}.rhs(1:startidx-1)) '*' ...
                getStrMoveRight(jsonmodel{i}.rhs(endidx+1:end))];
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
            error('Shouldn''t arrive here');
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
    if ~all(dropvname)
        vnames = vnames(dropvname);
    end

    lhssub = getRhsToSubFromLhs(ds, jsonmodel{i}.rhs, regex, splitstrings, pnames);

    residnames{i} = setdiff(intersect(rhs_, M_.exo_names), ds.name);
    assert(~isempty(residnames{i}), ['No residuals in equation ' num2str(i)]);
    assert(length(residnames{i}) == 1, ['More than one residual in equation ' num2str(i)]);

    params{i} = pnames;
    vars{i} = vnames;

    ydata = eval(regexprep(jsonmodel{i}.lhs, regex, 'ds.$&'));
    if ~isempty(lhssub)
        ydata = ydata - lhssub;
    end

    if isempty(xjdata)
        % AR(1) case
        fp = ydata.firstobservedperiod;
        lp = ydata.lastobservedperiod;
        startidxs(i) = length(Y) + 1;
        startdates{i} = fp;
        enddates{i} = lp;
        Y(startidxs(i):startidxs(i)+lp-fp, 1) = ydata(fp:lp).data;
        if columns(X) == 0
            X(startidxs(i):startidxs(i)+lp-fp, :) = zeros(ydata(fp:lp).nobs, 1);
        else
            X(startidxs(i):startidxs(i)+lp-fp, :) = zeros(ydata(fp:lp).nobs, columns(X));
        end
    else
        fp = max(ydata.firstobservedperiod, xjdata.firstobservedperiod);
        lp = min(ydata.lastobservedperiod, xjdata.lastobservedperiod);
        
        startidxs(i) = length(Y) + 1;
        startdates{i} = fp;
        enddates{i} = lp;
        Y(startidxs(i):startidxs(i)+lp-fp, 1) = ydata(fp:lp).data;
        X(startidxs(i):startidxs(i)+lp-fp, pidxs) = xjdata(fp:lp).data;
    end
end
surpidxs = surpidxs(1:surpidx);
end