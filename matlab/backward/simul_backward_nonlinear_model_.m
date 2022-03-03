function [ysim, xsim] = simul_backward_nonlinear_model_(initialconditions, samplesize, DynareOptions, DynareModel, DynareOutput, innovations, iy1, model_dynamic)

% Simulates a stochastic non linear backward looking model with arbitrary precision (a deterministic solver is used).
%
% INPUTS
% - initial_conditions  [dseries]     initial conditions for the endogenous variables.
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

% Copyright © 2017-2022 Dynare Team
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
        dprintf('Period t = %s.', num2str(it-initialconditions.nobs));
    end
    y_ = DynareOutput.endo_simul(:,it-1);
    ylag = y_(iy1);                    % Set lagged variables.
    y = y_;                            % A good guess for the initial conditions is the previous values for the endogenous variables.
    try
        if ismember(DynareOptions.solve_algo, [12,14])
            [DynareOutput.endo_simul(:,it), errorflag] = ...
                dynare_solve(model_dynamic_s, y, DynareOptions, ...
                             DynareModel.isloggedlhs, DynareModel.isauxdiffloggedrhs, DynareModel.endo_names, DynareModel.lhs, ...
                             model_dynamic, ylag, DynareOutput.exo_simul, DynareModel.params, DynareOutput.steady_state, it);
        else
            [DynareOutput.endo_simul(:,it), errorflag] = ...
                dynare_solve(model_dynamic_s, y, DynareOptions, ...
                             model_dynamic, ylag, DynareOutput.exo_simul, DynareModel.params, DynareOutput.steady_state, it);
        end
        if errorflag
            error()
        end
    catch
        DynareOutput.endo_simul = DynareOutput.endo_simul(:, 1:it-1);
        dprintf('Newton failed on iteration i = %s.', num2str(it-initialconditions.nobs));
        ytm = DynareOutput.endo_simul(:,end);
        xtt = DynareOutput.exo_simul(it,:);
        skipline()
        dprintf('Values of the endogenous variables before the nonlinear solver failure')
        dprintf('----------------------------------------------------------------------')
        skipline()
        dyntable(DynareOptions, '', {'VARIABLES','VALUES'}, DynareModel.endo_names(1:DynareModel.orig_endo_nbr), ytm(1:DynareModel.orig_endo_nbr), [], [], 6)
        skipline()
        dprintf('Values of the exogenous variables before the nonlinear solver failure')
        dprintf('---------------------------------------------------------------------')
        skipline()
        dyntable(DynareOptions, '', {'VARIABLES','VALUES'}, DynareModel.exo_names, transpose(DynareOutput.exo_simul(it,:)), [], [], 6)
        skipline(2)
        %
        % Get equation tags if any
        %
        if isfield(DynareModel, 'equations_tags')
            etags = cell(DynareModel.orig_endo_nbr, 1);
            for i = 1:DynareModel.orig_endo_nbr
                equations_tags = DynareModel.equations_tags(cellfun(@(x) isequal(x, i), DynareModel.equations_tags(:,1)), :);
                name = equations_tags(strcmpi(equations_tags(:,2), 'name'),:);
                if isempty(name)
                    eqtags{i} = int2str(i);
                else
                    if rows(name)>1
                        error('Something is wrong in the equation tags.')
                    else
                        eqtags(i) = name(3);
                    end
                end
            end
        else
            etags = split(int2str(1:DynareModel.orig_endo_nbr), '  ');
        end
        %
        % Evaluate and check the residuals
        %
        r = feval(model_dynamic, [ytm; ytm(iy1)], DynareOutput.exo_simul, DynareModel.params, DynareOutput.steady_state, it);
        residuals_evaluating_to_nan = isnan(r);
        residuals_evaluating_to_inf = isinf(r);
        residuals_evaluating_to_complex = ~isreal(r);
        if any(residuals_evaluating_to_nan)
            dprintf('Following equations are evaluating to NaN:')
            skipline()
            display_names_of_problematic_equations(DynareModel, residuals_evaluating_to_nan);
            skipline()
        end
        if any(residuals_evaluating_to_inf)
            dprintf('Following equations are evaluating to Inf:')
            skipline()
            display_names_of_problematic_equations(DynareModel, residuals_evaluating_to_nan);
            skipline()
        end
        if any(residuals_evaluating_to_complex)
            dprintf('Following equations are evaluating to a complex number:')
            skipline()
            display_names_of_problematic_equations(DynareModel, residuals_evaluating_to_nan);
            skipline()
        end
        % TODO Implement same checks with the jacobian matrix.
        break
    end
end

ysim = DynareOutput.endo_simul(1:DynareModel.orig_endo_nbr,:);
xsim = DynareOutput.exo_simul;


function display_names_of_problematic_equations(DynareModel, TruthTable)
for i=1:DynareModel.orig_endo_nbr
    if TruthTable(i)
        dprintf(' - %s', eqtags{i})
    end
end
for i=DynareModel.orig_endo_nbr+1:DynareModel.endo_nbr
    if TruthTable(i)
        dprintf(' - Auxiliary equation for %s', DynareModel.endo_names{i})
    end
end