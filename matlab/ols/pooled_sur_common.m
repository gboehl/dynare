function [X, Y, startdates, enddates, startidxs, residnames, pbeta, vars, surpidxs, surconstrainedparams] = pooled_sur_common(ds, lhs, rhs, lineno, M_exo_names_trim, M_param_names_trim);

%
% Code common to sur.m and pooled_ols.m
%
% INPUTS
%   none
%
% OUTPUTS
%   none
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

global M_

M_endo_exo_names_trim = [cellstr(M_.endo_names); M_exo_names_trim];
regex = strjoin(M_endo_exo_names_trim(:,1), '|');
mathops = '[\+\*\^\-\/]';
params = cell(length(rhs),1);
vars = cell(length(rhs),1);
pbeta = {};
Y = [];
X = [];
startidxs = zeros(length(lhs), 1);
startdates = cell(length(lhs), 1);
enddates = cell(length(lhs), 1);
residnames = cell(length(lhs), 1);
surpidxs = zeros(M_.param_nbr, 1); 
surpidx = 0; 
surconstrainedparams = [];
for i = 1:length(lhs)
    rhs_ = strsplit(rhs{i}, {'+','-','*','/','^','log(','ln(','log10(','exp(','(',')','diff('});
    rhs_(cellfun(@(x) all(isstrprop(x, 'digit')), rhs_)) = [];
    vnames = setdiff(rhs_, M_param_names_trim);
    if ~isempty(regexp(rhs{i}, ...
            ['(' strjoin(vnames, '\\(\\d+\\)|') '\\(\\d+\\))'], 'once'))
        error(['pooled_ols: you cannot have leads in equation on line ' ...
            lineno{i} ': ' lhs{i} ' = ' rhs{i}]);
    end

    % Find parameters and associated variables
    pnames = intersect(rhs_, M_param_names_trim);
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
            surpidxs(surpidx, 1) = find(strcmp(pnames{j}, M_param_names_trim));
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
        [startidx, endidx] = regexp(rhs{i}, pregex, 'start', 'end');
        assert(length(startidx) == 1);
        if rhs{i}(startidx) == '*' && rhs{i}(endidx) == '*'
            vnames{j} = [getStrMoveLeft(rhs{i}(1:startidx-1)) '*' ...
                getStrMoveRight(rhs{i}(endidx+1:end))];
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
            error('pooled_ols: Shouldn''t arrive here');
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

    lhssub = getRhsToSubFromLhs(ds, rhs{i}, regex, [splitstrings; pnames]);

    residnames{i} = setdiff(intersect(rhs_, M_exo_names_trim), ds.name);
    assert(~isempty(residnames{i}), ['No residuals in equation ' num2str(i)]);
    assert(length(residnames{i}) == 1, ['More than one residual in equation ' num2str(i)]);

    params{i} = pnames;
    vars{i} = vnames;

    ydata = eval(regexprep(lhs{i}, regex, 'ds.$&'));
    for j = 1:lhssub.vobs
        ydata = ydata - lhssub{j};
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