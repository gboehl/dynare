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

% Set up initial conditions
[initialcondition, periods, innovations, DynareOptions, DynareModel, DynareOutput, endonames, exonames, nx, ny1, iy1, jdx, model_dynamic, y] = ...
    simul_backward_model_init(initialcondition, periods, options_, M_, oo_, zeros(periods, M_.exo_nbr));

% Get the covariance matrix of the shocks.
Sigma = M_.Sigma_e + 1e-14*eye(M_.exo_nbr);
sigma = transpose(chol(Sigma));

% Initialization of the returned argument. Each will be a dseries object containing the IRFS for the endogenous variables listed in the third input argument.
irfs = struct();

% Compute the IRFs (loop over innovations).
for i=1:length(listofshocks)
    innovations = zeros(periods, DynareModel.exo_nbr);
    % Get transition paths induced by the initial condition.
    if options_.linear
        ysim__0 = simul_backward_linear_model_(initialcondition, periods, DynareOptions, DynareModel, DynareOutput, innovations, nx, ny1, iy1, jdx, model_dynamic);
    else
        ysim__0 = simul_backward_nonlinear_model_(initialcondition, periods, DynareOptions, DynareModel, DynareOutput, innovations, iy1, model_dynamic);
    end
    % Add the shock.
    j = find(strcmp(listofshocks{i}, exonames));
    if isempty(j)
        error('backward_model_irf: Exogenous variable %s is unknown!', listofshocks{i})
    end
    innovations(1,:) = transpose(sigma(:,j));
    if options_.linear
        ysim__1 = simul_backward_linear_model_(initialcondition, periods, DynareOptions, DynareModel, DynareOutput, innovations, nx, ny1, iy1, jdx, model_dynamic);
    else
        ysim__1 = simul_backward_nonlinear_model_(initialcondition, periods, DynareOptions, DynareModel, DynareOutput, innovations, iy1, model_dynamic);
    end
    % Transform the endogenous variables
    if notransform
        endo_simul__0 = ysim__0;
        endo_simul__1 = ysim__1;
    else
        endo_simul__0 = feval(transform, ysim__0);
        endo_simul__1 = feval(transform, ysim__1);
    end
    % Instantiate a dseries object (with all the endogenous variables)
    allirfs = dseries(transpose(endo_simul__1-endo_simul__0), initialcondition.init, endonames, cellstr(DynareModel.endo_names_tex));
    % Extract a sub-dseries object
    irfs.(listofshocks{i}) = allirfs{listofvariables{:}};
end