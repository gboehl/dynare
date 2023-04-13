function [ysim, xsim, errorflag] = simul_backward_linear_model_(initialconditions, samplesize, DynareOptions, DynareModel, DynareOutput, innovations, dynamic_resid, dynamic_g1)

% Simulates a stochastic linear backward looking model.
%
% INPUTS
% - initialconditions   [dseries]     initial conditions for the endogenous variables.
% - samplesize          [integer]     scalar, number of periods for the simulation.
% - DynareOptions       [struct]      Dynare's options_ global structure.
% - DynareModel         [struct]      Dynare's M_ global structure.
% - DynareOutput        [struct]      Dynare's oo_ global structure.
% - innovations         [double]      T*q matrix, innovations to be used for the simulation.
%
% OUTPUTS
% - DynareOutput        [struct]      Dynare's oo_ global structure.
% - errorflag           [logical]     scalar, equal to false iff the simulation did not fail.
%
% REMARKS
% [1] The innovations used for the simulation are saved in DynareOutput.exo_simul, and the resulting paths for the endogenous
%     variables are saved in DynareOutput.endo_simul.
% [2] The last input argument is not mandatory. If absent we use random draws and rescale them with the informations provided
%     through the shocks block.
% [3] If the first input argument is empty, the endogenous variables are initialized with 0, or if available with the informations
%     provided thrtough the histval block.

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

errorflag = false;

if ~isempty(innovations)
    DynareOutput.exo_simul(initialconditions.nobs+(1:samplesize),:) = innovations;
end

% Get coefficients
y = [zeros(2*DynareModel.endo_nbr,1); NaN(DynareModel.endo_nbr,1)];
x = zeros(DynareModel.exo_nbr, 1);
[cst, T_order, T] = dynamic_resid(y, x, DynareModel.params, DynareOutput.steady_state);
jacob = dynamic_g1(y, x, DynareModel.params, DynareOutput.steady_state, DynareModel.dynamic_g1_sparse_rowval, DynareModel.dynamic_g1_sparse_colval, DynareModel.dynamic_g1_sparse_colptr, T_order, T);

try
    A0inv = inv(jacob(:,DynareModel.endo_nbr+(1:DynareModel.endo_nbr)));
catch
    errorflag = true;
    ysim = [];
    xsim = [];
    return
end

A1 = jacob(:,1:DynareModel.endo_nbr);
B = jacob(:,3*DynareModel.endo_nbr+1:end);

% Simulations
for it = initialconditions.nobs+(1:samplesize)
    DynareOutput.endo_simul(:,it) = -A0inv*(cst + A1*DynareOutput.endo_simul(:,it-1) + B*DynareOutput.exo_simul(it,:)');
end

ysim = DynareOutput.endo_simul(1:DynareModel.orig_endo_nbr,:);
xsim = DynareOutput.exo_simul;
