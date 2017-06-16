function irfs = backward_model_irf(initialcondition, listofshocks, listofvariables, varargin)

% Returns impulse response functions.
%
% INPUTS 
% - initialcondition    [dseries,dates]       Initial conditions for the endogenous variables, or period 0.
% - listofshocks        [cell of strings]     The innovations for which the IRFs need to be computed.  
% - listofvariables     [cell of strings]     The endogenous variables which will be returned.  
% - periods             [integer]             scalar, the number of periods.
%
% OUTPUTS 
% - irfs                [struct of dseries]   
%
% REMARKS 
% The names of the fields in the returned structure are given by the name
% of the innovations listed in the second input argument. Each field gather
% the associated paths for endogenous variables listed in the third input 
% argument. 


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
    error(['simul_model_irf:: The specified model is not backward looking!'])
end

% Set default value for the fourth input argument.
if nargin<4
    periods = 40;
    notransform = true;
else
    periods = varargin{1};
end

% Set default value for the last input argument (no transformation).
if nargin<5
    notransform = true;
else
    notransform = false;
    transform = varargin{2};
end

% Get list of all exogenous variables in a cell of strings
exo_names = cellstr(M_.exo_names);

% Get the list of all exogenous variables in a cell of strings
endo_names = cellstr(M_.endo_names);

% Set initial condition.
if isdates(initialcondition)
    if isempty(M_.endo_histval)
        error('backward_model_irf: histval block for setting initial condition is missing!')
    end
    initialcondition = dseries(transpose(M_.endo_histval), initialcondition, endo_names, cellstr(M_.endo_names_tex));
end

% Get the covariance matrix of the shocks.
Sigma = M_.Sigma_e + 1e-14*eye(M_.exo_nbr);
sigma = transpose(chol(Sigma));

% Put initial conditions in a vector of doubles
initialconditions = transpose(initialcondition{endo_names{:}}.data);

% Initialization of the returned argument. Each will be a dseries object containing the IRFS for the endogenous variables listed in the third input argument.
irfs = struct();

% Get the covariance matrix of the shocks.
Sigma = M_.Sigma_e + 1e-14*eye(M_.exo_nbr);
sigma = transpose(chol(Sigma));

% Put initial conditions in a vector of doubles
initialconditions = transpose(initialcondition{endo_names{:}}.data);

% Compute the IRFs (loop over innovations).
for i=1:length(listofshocks)
    % Get transition paths induced by the initial condition.
    innovations = zeros(periods+M_.maximum_exo_lag, M_.exo_nbr);
    if ~isempty(M_.exo_histval)
        innovations(1:M_.maximum_exo_lag,:) = M_.exo_histval;
    end
    oo__0 = simul_backward_model(initialconditions, periods, options_, M_, oo_, innovations);
    % Add the shock.
    j = strmatch(listofshocks{i}, exo_names);
    if isempty(j)
        error('backward_model_irf: Exogenous variable %s is unknown!', listofshocks{i})
    end
    innovations(1+M_.maximum_exo_lag,:) = transpose(sigma(:,j));
    oo__1 = simul_backward_model(initialconditions, periods, options_, M_, oo_, innovations);
    % Transform the endogenous variables
    if notransform
        endo_simul__0 = oo__0.endo_simul;
        endo_simul__1 = oo__1.endo_simul;
    else
        endo_simul__0 = feval(transform, oo__0.endo_simul);
        endo_simul__1 = feval(transform, oo__1.endo_simul);
    end
    % Instantiate a dseries object (with all the endogenous variables)
    allirfs = dseries(transpose(endo_simul__1-endo_simul__0), initialcondition.init, cellstr(M_.endo_names), cellstr(M_.endo_names_tex));
    % Extract a sub-dseries object
    irfs.(listofshocks{i}) = allirfs{listofvariables{:}};
end