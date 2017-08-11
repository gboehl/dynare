function forecasts = backward_model_forecast(initialcondition, listofvariables, periods, withuncertainty)

% Returns unconditional forecasts.
%
% INPUTS 
% - initialcondition    [dseries]             Initial conditions for the endogenous variables.
% - periods             [integer]             scalar, the number of (forecast) periods.
% - withuncertainty     [logical]             scalar, returns confidence bands if true.
%
% OUTPUTS 
% - forecast            [dseries]   

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

global M_ options_ oo_

% Check that the model is actually backward
if M_.maximum_lead
    error(['backward_model_irf:: The specified model is not backward looking!'])
end

% Initialize returned argument.
forecasts = struct();

% Set defaults.
if nargin<2
    listofvariables = cellstr(M_.endo_names);
    periods = 8;
    withuncertainty = false;
end

if nargin<3
    periods = 8;
    withuncertainty = false;
end

if nargin<4
    withuncertainty = false;
end

% Get full list of endogenous variables
endo_names = cellstr(M_.endo_names);

% Get vector of indices for the selected endogenous variables.
n = length(listofvariables);
idy = zeros(n,1);
for i=1:n
    j = find(strcmp(listofvariables{i}, endo_names));
    if isempty(j)
        error('backward_model_forecast:: Variable %s is unknown!', listofvariables{i})
    else
        idy(i) = j;
    end
end

% Set the number of simulations (if required).
if withuncertainty
    B = 1000;
end

% Get the covariance matrix of the shocks.
if withuncertainty
    Sigma = M_.Sigma_e + 1e-14*eye(M_.exo_nbr);
    sigma = transpose(chol(Sigma));
end

% Set initial condition.
if isdates(initialcondition)
    if isempty(M_.endo_histval)
        error('backward_model_irf: histval block for setting initial condition is missing!')
    end
    initialcondition = dseries(transpose(M_.endo_histval), initialcondition, endo_names, cellstr(M_.endo_names_tex));
end

% Put initial conditions in a vector of doubles
initialconditions = transpose(initialcondition{endo_names{:}}.data);

% Compute forecast without shock 
innovations = zeros(periods+M_.maximum_exo_lag, M_.exo_nbr);
if M_.maximum_exo_lag
    if isempty(M_.exo_histval)
        error('You need to set the past values for the exogenous variables!')
    else
        innovations(1:M_.maximum_exo_lag, :) = M_.exo_histval;
    end
end

oo__0 = simul_backward_model(initialconditions, periods, options_, M_, oo_, innovations);
forecasts.pointforecast = dseries(transpose(oo__0.endo_simul(idy,:)), initialcondition.init, listofvariables);

if withuncertainty
    % Preallocate an array gathering the simulations.
    ArrayOfForecasts = zeros(n, periods+size(initialconditions, 2), B);
    for i=1:B
        innovations(M_.maximum_exo_lag+1:end,:) = transpose(sigma*randn(M_.exo_nbr, periods));
        oo__ = simul_backward_model(initialconditions, periods, options_, M_, oo_, innovations);     
        ArrayOfForecasts(:,:,i) = oo__.endo_simul(idy,:); 
    end
    % Compute mean (over future uncertainty) forecast.
    forecasts.meanforecast = dseries(transpose(mean(ArrayOfForecasts, 3)), initialcondition.init, listofvariables);
    forecasts.medianforecast = dseries(transpose(median(ArrayOfForecasts, 3)), initialcondition.init, listofvariables);
    forecasts.stdforecast = dseries(transpose(std(ArrayOfForecasts, 1,3)), initialcondition.init, listofvariables);
    % Compute lower and upper 95% confidence bands
    ArrayOfForecasts = sort(ArrayOfForecasts, 3);
    forecasts.lb = dseries(transpose(ArrayOfForecasts(:,:,round(0.025*B))), initialcondition.init, listofvariables);
    forecasts.ub = dseries(transpose(ArrayOfForecasts(:,:,round(0.975*B))), initialcondition.init, listofvariables);
end