function plot_contributions(tagn, tagv, dseriesdata, params)
% function plot_contributions(tagn, tagv, dseriesdata, params)
% Return the lhs, rhs of an equation and the line it was defined
% on given its tag
%
% INPUTS
%   tagn           string        tag name
%   tagv           string        tag value
%   dseriesdata    matrix        nobs x n endogenous in equation
%   params         struct        params.param_name = param_value
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

global M_;

jsonfile = [M_.fname '_original.json'];
if exist(jsonfile, 'file') ~= 2
    error(['Could not find ' jsonfile]);
end

% Get equation
jsonmodel = loadjson(jsonfile);
jsonmodel = jsonmodel.model;
[lhs, rhs, ~] = getEquationByTag(jsonmodel, tagn, tagv);

w% replace variables with dseriesdata.variablename
for i = 1:length(dseriesdata.name)
    rhs = strrep(rhs, [dseriesdata{i}.name{:} '('], ['dseriesdata.' dseriesdata{i}.name{:} '(']);
end

% create inline function for rhs
rhsfunc = inline(rhs);

% construct call
rhsfuncargs = argnames(rhsfunc);
nargs = length(rhsfuncargs);
parameters = cell(nargs-1, 1);
for i = 2:nargs
    parameters{i-1} = params.(rhsfuncargs{i});
end

contribution = zeros(dseriesdata.nobs, nargs+1);
funceval = rhsfunc(dseriesdata, parameters{:});
contribution(:, 1) = funceval.data;

nvars = length(dseriesdata.name);
dseriesdatabak = dseriesdata;
dseriesdatazero = dseries(zeros(dseriesdata.nobs, nvars), ...
    dseriesdata.firstdate, ...
    dseriesdata.name);
for i = 1:nvars
    dseriesdata = dseriesdatazero;
    dseriesdata{dseriesdatabak.name{i}} = dseriesdatabak{i};
    funceval = rhsfunc(dseriesdata, parameters{:});
    contribution(:, i+1) = funceval.data;
end

figure('Name', lhs);
plot(1:dseriesdata.nobs, contribution)
seriesnames = dseriesdata.name;
legend('All Endogenous',seriesnames{:})

end