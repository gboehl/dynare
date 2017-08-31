function DynareOutput = simul_backward_nonlinear_model(varargin)

% Simulates a stochastic non linear backward looking model with arbitrary precision (a deterministic solver is used).
%
% INPUTS
% - initial_conditions  [double]      n*1 vector, initial conditions for the endogenous variables.
% - sample_size         [integer]     scalar, number of periods for the simulation.
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

model_dynamic_s = str2func('dynamic_backward_model_for_simulation');

% initialization of the returned simulations.
DynareOutput.endo_simul = NaN(DynareModel.endo_nbr, samplesize);
if isempty(initialconditions)
    if isfield(DynareModel,'endo_histval') && ~isempty(DynareModel.endo_histval)
        DynareOutput.endo_simul = [DynareModel.endo_histval, DynareOutput.endo_simul];
    else
        warning('simul_backward_nonlinear_model:: Initial condition is zero for all variables! The model simulation may fail with the default initialization.')
        DynareOutput.endo_simul = [zeros(DynareModel.endo_nbr, DynareModel.max_lag_orig), DynareOutput.endo_simul];
    end
else
    if ~isequal(size(initialconditions, 2), DynareModel.max_lag_orig)
        error(['simul_backward_nonlinear_model:: First argument should have %s columns!'], DynareModel.max_lag_orig)
    end
    DynareOutput.endo_simul = [initialconditions, DynareOutput.endo_simul];
end
Y = DynareOutput.endo_simul;

DynareOutput.exo_simul = [ zeros(1, DynareModel.exo_nbr); DynareOutput.exo_simul];

% Simulations (call a Newton-like algorithm for each period).
for it = 1+(1:samplesize)
    ylag = Y(iy1,it-1);                   % Set lagged variables.
    y = Y(:,it-1);                        % A good guess for the initial conditions is the previous values for the endogenous variables.
    Y(:,it) = dynare_solve(model_dynamic_s, y, DynareOptions, model_dynamic, ylag, DynareOutput.exo_simul, DynareModel.params, DynareOutput.steady_state, it);
end

DynareOutput.endo_simul = Y;