function [residuals, info] = calibrateresiduals(dbase, info, DynareModel)

% Compute residuals in a backward model. Residuals are unobserved exogenous
% variables appearing additively in equations and without lags. An equation
% cannot have more than one residual, and a residual cannot appear in more
% than one equation.
%
% INPUTS
% - dbase       [dseries]   Object containing all the endogenous and observed exogenous variables.
% - info        [struct]    Informations about the residuals.
% - DynareModel [struct]    M_ as produced by the preprocessor.
%
% OUTPUTS
% - residuals   [dseries]   Object containing the identified residuals.
% - info        [struct]    Informations about the residuals.
%
% REMARKS
% The first two input arguments are the output of checkdatabaseforinversion
% routine.

% Copyright © 2017-2020 Dynare Team
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
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

displayresidualsequationmapping = false;

% Get function handle for the dynamic model
model_dynamic = str2func([DynareModel.fname,'.dynamic']);

% Get data for all the endogenous variables.
ydata = dbase{info.endonames{:}}.data;

% Define function to retrieve an equation name
eqname = @(z) DynareModel.equations_tags{cellfun(@(x) x==z, DynareModel.equations_tags(:,1)) & cellfun(@(x) isequal(x, 'name'), DynareModel.equations_tags(:,2)),3};

% Get data for all the exogenous variables. Missing exogenous variables, to be solved for, have NaN values.
exogenousvariablesindbase = intersect(info.exonames, dbase.name);
residuals = dseries(NaN(dbase.nobs, length(info.residuals)), dbase.init, info.residuals);
allexogenousvariables = [dbase{exogenousvariablesindbase{:}}, residuals];
allexogenousvariables = allexogenousvariables{info.exonames{:}};
xdata = allexogenousvariables.data;

% Evaluate the dynamic equation
n = size(ydata, 2);
c = find(DynareModel.lead_lag_incidence');
y = [ydata(1,:)'; ydata(2,:)'];
y = y(c);
r = model_dynamic(y, xdata, DynareModel.params, zeros(n, 1), 2);

% Check that the number of equations evaluating to NaN matches the number of residuals
idr = find(isnan(r));
if ~isequal(length(idr), residuals.vobs)
    error('Each residual should appear in only one equation, and an equation cannot have more than one residual!')
end

% Check that the non NaN equations have zero residuals (model and data consistency).
ido = setdiff(1:n, idr);
if ~isempty(find(abs(r(ido))>1e-6))
    disp('Provided data and model are not consistent in equations:')
    idx = find(abs(r)>1e-6);
    c1 = 'Equation';
    c1 = strvcat(c1, '--------');
    for i = 1:length(idx)
        c1 = strvcat(c1, sprintf('  %s', num2str(idx(i))));
    end
    c2 = 'Residual';
    c2 = strvcat(c2, '--------');
    for i = 1:length(idx)
        c2 = strvcat(c2, sprintf('%s', num2str(r(idx(i)))));
    end
    c2 = strvcat(c2(1, :), repmat('-', 1, size(c2, 2)), c2(3:end,:));
    c5 = 'Equation name';
    c5 = strvcat(c5, '-------------');
    for i = 1:length(idx)
        c5 = strvcat(c5, sprintf('  %s', eqname(idx(i))));
    end
    c3 = repmat(' | ', size(c2, 1), 1);
    c4 = repmat('   ', size(c2, 1), 1);
    cc = [c4, c1, c3, c2, c3, c5];
    skipline()
    disp(cc)
    skipline()
    disp('Please check model and dataset.')
end

% Associate the residuals with equations equations evaluating to NaNs.
info.equations = cell(residuals.vobs, 1);
info.residualindex = cell(residuals.vobs, 1);
for i = 1:residuals.vobs
    residualname = residuals.name{i};
    info.residualindex(i) = {strmatch(residualname, allexogenousvariables.name, 'exact')};
    tmpxdata = xdata;
    tmpxdata(2, info.residualindex{i}) = 0;
    r = model_dynamic(y, tmpxdata, DynareModel.params, zeros(n, 1), 2);
    info.equations(i) = { idr(find(~isnan(r(idr))))};
end

if displayresidualsequationmapping
    c1 = 'Residual';
    for i=1:length(info.residuals)
        c1 = strvcat(c1, sprintf('%s', info.residuals{i}));
    end
    c1 = strvcat(c1(1,:), repmat('-', 1, size(c1, 2)), c1(2:end,:));
    c2 = 'Equation';
    for i=1:length(info.residuals)
        c2 = strvcat(c2, sprintf('  %s', num2str(info.equations{i})));
    end
    c2 = strvcat(c2(1,:), repmat('-', 1, size(c2, 2)), c2(2:end,:));
    c3 = repmat(' | ', size(c2, 1), 1);
    c4 = repmat('   ', size(c2, 1), 1);
    cc = [c4, c1, c3, c2];
    skipline()
    disp(cc)
    skipline()
end

% Compute residuals
xdata(:,cell2mat(info.residualindex)) = 0;
rdata = NaN(residuals.nobs, residuals.vobs);
for t=2:size(xdata, 1)
    y = transpose([ydata(t-1,:); ydata(t,:)]);
    r = model_dynamic(y(c), xdata, DynareModel.params, zeros(n, 1), t);
    rdata(t,:) = transpose(r(cell2mat(info.equations)));
end
residuals = dseries(rdata, dbase.init, info.residuals);