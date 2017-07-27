function DynareOutput = simul_backward_linear_model(varargin)

% Simulates a stochastic linear backward looking model.
%
% INPUTS
% - initialconditions   [double]      n*1 vector, initial conditions for the endogenous variables.
% - samplesize          [integer]     scalar, number of periods for the simulation.
% - DynareOptions       [struct]      Dynare's options_ global structure.
% - DynareModel         [struct]      Dynare's M_ global structure.
% - DynareOutput        [struct]      Dynare's oo_ global structure.
% - innovations         [double]      T*q matrix, innovations to be used for the simulation.
%
% OUTPUTS
% - DynareOutput        [struct]      Dynare's oo_ global structure.
%
% REMARKS
% [1] The innovations used for the simulation are saved in DynareOutput.exo_simul, and the resulting paths for the endogenous
%     variables are saved in DynareOutput.endo_simul.
% [2] The last input argument is not mandatory. If absent we use random draws and rescale them with the informations provided
%     through the shocks block.
% [3] If the first input argument is empty, the endogenous variables are initialized with 0, or if available with the informations
%     provided thrtough the histval block.

% Copyright (C) 2012-2017 Dynare Team
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

[initialconditions, samplesize, innovations, DynareOptions, DynareModel, DynareOutput, nx, ny1, iy1, jdx, model_dynamic, y] = ...
    simul_backward_model_init(varargin{:});

% initialization of the returned simulations.
DynareOutput.endo_simul = NaN(DynareModel.endo_nbr,samplesize+1);
if isempty(initialconditions)
    DynareOutput.endo_simul(:,1) = DynareOutput.steady_state;
else
    DynareOutput.endo_simul(:,1) = initialconditions;
end
Y = DynareOutput.endo_simul;

% get coefficients
[cst,jacob] = model_dynamic(zeros(DynareModel.endo_nbr+ny1,1), ...
                            zeros(2,DynareModel.exo_nbr), ...
                            DynareModel.params, ...
                            DynareOutput.steady_state,1);
A0inv = inv(jacob(:,jdx));
A1 = jacob(:,nonzeros(DynareModel.lead_lag_incidence(1,:)));
B = jacob(:,end-nx+1:end);

% Simulations
for it = 2:samplesize+1
    Y(:,it) = -A0inv*(cst + A1*Y(iy1,it-1) + B*DynareOutput.exo_simul(it,:)');
end

DynareOutput.endo_simul = Y;