function forecasts = backward_model_forecast(initialcondition, listofvariables, periods, withuncertainty)
%function forecasts = backward_model_forecast(initialcondition, listofvariables, periods, withuncertainty)
% Returns unconditional forecasts.
%
% INPUTS
% - initialcondition    [dseries]             Initial conditions for the endogenous variables.
% - periods             [integer]             scalar, the number of (forecast) periods.
% - withuncertainty     [logical]             scalar, returns confidence bands if true.
%
% OUTPUTS
% - forecast            [dseries]

% Copyright © 2017-2023 Dynare Team
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

global M_ options_ oo_

% Check that the model is actually backward
if M_.maximum_lead
    error('backward_model_forecast:: The specified model is not backward looking!')
end

% Initialize returned argument.
forecasts = struct();

% Set defaults.
if nargin<2
    listofvariables = M_.endo_names;
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

start = initialcondition.dates(end)+1;

% Set default initial conditions for the innovations.
for i=1:M_.exo_nbr
    if ~ismember(M_.exo_names{i}, initialcondition.name)
        initialcondition{M_.exo_names{i}} = dseries(zeros(initialcondition.nobs, 1), initialcondition.dates(1), M_.exo_names{i});
    end
end

% Set up initial conditions
[initialcondition, periods, innovations, options_local, M_local, oo_local, endonames, ~, dynamic_resid, dynamic_g1] = ...
    simul_backward_model_init(initialcondition, periods, options_, M_, oo_, zeros(periods, M_.exo_nbr));

% Get vector of indices for the selected endogenous variables.
n = length(listofvariables);
idy = zeros(n,1);
for i=1:n
    j = find(strcmp(listofvariables{i}, endonames));
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
    sigma = get_lower_cholesky_covariance(M_.Sigma_e,options_.add_tiny_number_to_cholesky);
end

% Compute forecast without shock
if options_.linear
    [ysim__0, errorflag] = simul_backward_linear_model_(initialcondition, periods, options_local, M_local, oo_local, innovations, dynamic_resid, dynamic_g1);
else
    [ysim__0, errorflag] = simul_backward_nonlinear_model_(initialcondition, periods, options_local, M_local, oo_local, innovations, dynamic_resid, dynamic_g1);
end

if errorflag
    error('Simulation failed.')
end

forecasts.pointforecast = dseries(transpose(ysim__0(idy,:)), initialcondition.init, listofvariables);

% Set first period of forecast
forecasts.start = start;

if withuncertainty
    % Preallocate an array gathering the simulations.
    ArrayOfForecasts = zeros(n, periods+initialcondition.nobs, B);
    for i=1:B
        innovations = transpose(sigma*randn(M_.exo_nbr, periods));
        if options_.linear
            [ysim__, ~, errorflag] = simul_backward_linear_model_(initialcondition, periods, options_local, M_local, oo_local, innovations, dynamic_resid, dynamic_g1);
        else
            [ysim__, ~, errorflag] = simul_backward_nonlinear_model_(initialcondition, periods, options_local, M_local, oo_local, innovations, dynamic_resid, dynamic_g1);
        end
        if errorflag
            error('Simulation failed.')
        end
        ArrayOfForecasts(:,:,i) = ysim__(idy,:);
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
