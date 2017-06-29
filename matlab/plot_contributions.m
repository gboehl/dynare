function plot_contributions(tagn, tagv, ds1, ds0)

% Plots the contribution to the lhs variable of the rhs variables in an equation.
%
% INPUTS
%  - tagn      [string]                 Equation tag name
%  - tagv      [string]                 Equation tag value
%  - ds1       [dseries]                Object containing all the variables (exogenous and endogenous) appearing in the equation.
%  - ds0       [dseries]                parameter values
%
% OUTPUTS
%   None
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

jsonfile = [M_.fname '_original.json'];
if exist(jsonfile, 'file') ~= 2
    error('Could not find %s! Please use the json option (See the Dynare invocation section in the reference manual).', jsonfile);
end

% Get equation.
jsonmodel = loadjson(jsonfile);
jsonmodel = jsonmodel.model;
[lhs, rhs, ~] = getEquationByTag(jsonmodel, tagn, tagv);

% Get variable and parameter names in the equation.
rhs_ = strsplit(rhs,{'+','-','*','/','^','log(','exp(','(',')'});
rhs_(find(cellfun(@(x) all(isstrprop(x, 'digit')), rhs_))) = []; % Remove numbers
pnames = cellstr(M_.param_names);
vnames = setdiff(rhs_, pnames);
pnames = setdiff(rhs_, vnames);

% Get values for the parameters
idp = strmatch(pnames{1}, M_.param_names, 'exact');
str = sprintf('%s = M_.params(%d);', pnames{1}, idp);
for i=2:length(pnames)
    idp = strmatch(pnames{i}, M_.param_names, 'exact');
    str = sprintf('%s %s = M_.params(%d);', str, pnames{i}, idp);
end
eval(str)

rhs

% Replace variables with ds.variablename
for i = 1:length(vnames)
    if ismember(vnames{i}, ds1.name) && ismember(vnames{i}, ds0.name)
        regularexpression = ['(([\s\+\-\*\/\^\)])|^)(' vnames{i} ')(([\s\(\+\-\*\/\^])|$)'];
        newstring = ['$1ds.' vnames{i} '$3'];
        rhs = regexprep(rhs, regularexpression, newstring);
    else
        if ismember(vnames{i}, ds1.name)
            error('Variable %s is not available in the second dseries (baseline paths)!', vnames{i})
        else
            error('Variable %s is not available in the first dseries (actual paths)!', vnames{i})
        end
    end
end

%  Evaluate RHS with actual paths for the endogenous variables (this should
%  be equal to LHS member of the equation).
ds = ds1;
contribution = zeros(ds.nobs, ds.vobs + 1);
rhseval = eval(rhs);
contribution(:, 1) = rhseval.data;

% Compute the marginal effect of each variable on the RHS, by evaluating the
% RHS with all variables at the Baseline paths except one for which the
% actual path is used.
for i = 1:length(vnames)
    ds = ds0; % Set all variable to Baseline paths.
    ds{vnames{i}} = ds1{vnames{i}};
    rhseval = eval(rhs);
    contribution(:, i+1) = rhseval.data;
end

% Create the contributions plot.
figure('Name', lhs);
plot(1:ds.nobs, contribution);
seriesnames = ds.name;
legend('All Endogenous', seriesnames{:});