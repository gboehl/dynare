function plot_contributions(tagn, tagv, ds, params)

% Plots the contribution to the lhs variable of the rhs variables in an equation.
%
% INPUTS
%  - tagn      [string]                 Equation tag name
%  - tagv      [string]                 Equation tag value
%  - ds        [dseries]                Object containing all the variables (exogenous and endogenous) appearing in the equation.
%  - params    [struct]                 parameter values
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

% Get equation
jsonmodel = loadjson(jsonfile);
jsonmodel = jsonmodel.model;
[lhs, rhs, ~] = getEquationByTag(jsonmodel, tagn, tagv);

% replace variables with ds.variablename
for i = 1:length(ds.name)
    rhs = regexprep(rhs, ['([\s\+\-\*\/\^]{1}|^)(' ds{i}.name{:} ')([\s\(\+\-\*\/\^|$]{1})'], '$1ds.$2$3');
end
fields = fieldnames(params);

for i = 1:length(fields)
    rhs = regexprep(rhs, ['([\s\+\-\*\/\^]{1}|^)(' fields{i} ')([\s\+\-\*\/\^]{1}|$)'], '$1params.$2$3');
end

% call function with all variable values
contribution = zeros(ds.nobs, ds.vobs + 1);
rhseval = eval(rhs);
contribution(:, 1) = rhseval.data;

dsbak = ds;
dszero = dseries(zeros(ds.nobs, ds.vobs), ...
                          ds.firstdate, ...
                          ds.name);

for i = 1:ds.vobs
    ds = dszero;
    ds{dsbak.name{i}} = dsbak{dsbak.name{i}};
    rhseval = eval(rhs);
    contribution(:, i+1) = rhseval.data;
end

figure('Name', lhs);
plot(1:ds.nobs, contribution);
seriesnames = ds.name;
legend('All Endogenous', seriesnames{:});

end