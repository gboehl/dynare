function [ysim, xsim] = simul_backward_nonlinear_model_(initialconditions, samplesize, DynareOptions, DynareModel, DynareOutput, innovations, iy1, model_dynamic)

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

debug = false;

model_dynamic_s = str2func('dynamic_backward_model_for_simulation');

if ~isempty(innovations)
    DynareOutput.exo_simul(initialconditions.nobs+(1:samplesize),:) = innovations;
end

% Simulations (call a Newton-like algorithm for each period).
for it = initialconditions.nobs+(1:samplesize)
    if debug
        dprintf('Période t = %s.', num2str(it-initialconditions.nobs));
    end
    y_ = DynareOutput.endo_simul(:,it-1);
    ylag = y_(iy1);                    % Set lagged variables.
    y = y_;                            % A good guess for the initial conditions is the previous values for the endogenous variables.
    try
        if ismember(DynareOptions.solve_algo, [12,14])
            [DynareOutput.endo_simul(:,it), info] = ...
                dynare_solve(model_dynamic_s, y, DynareOptions, ...
                             DynareModel.isloggedlhs, DynareModel.isauxdiffloggedrhs, DynareModel.endo_names, DynareModel.lhs, ...
                             model_dynamic, ylag, DynareOutput.exo_simul, DynareModel.params, DynareOutput.steady_state, it);
        else
            [DynareOutput.endo_simul(:,it), info] = ...
                dynare_solve(model_dynamic_s, y, DynareOptions, ...
                             model_dynamic, ylag, DynareOutput.exo_simul, DynareModel.params, DynareOutput.steady_state, it);
        end
        if info
            error('Newton failed!')
        end
    catch
        DynareOutput.endo_simul(:, 1:it-1);
        dprintf('Newton failed on iteration i = %s.', num2str(it-initialconditions.nobs));
        break
    end
end

ysim = DynareOutput.endo_simul(1:DynareModel.orig_endo_nbr,:);
xsim = DynareOutput.exo_simul;