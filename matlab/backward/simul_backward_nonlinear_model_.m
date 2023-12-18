function [ysim, xsim, errorflag] = simul_backward_nonlinear_model_(initialconditions, samplesize, options_, M_, oo_, innovations, dynamic_resid, dynamic_g1)
% [ysim, xsim, errorflag] = simul_backward_nonlinear_model_(initialconditions, samplesize, options_, M_, oo_, innovations, dynamic_resid, dynamic_g1)
% Simulates a stochastic non linear backward looking model with arbitrary precision (a deterministic solver is used).
%
% INPUTS
% - initial_conditions  [dseries]     initial conditions for the endogenous variables.
% - sample_size         [integer]     scalar, number of periods for the simulation.
% - options_            [struct]      Dynare's options_ global structure.
% - M_                  [struct]      Dynare's M_ global structure.
% - oo_                 [struct]      Dynare's oo_ global structure.
% - innovations         [double]      T*q matrix, innovations to be used for the simulation.
% - dynamic_resid
% - dynamic_g1
%
% OUTPUTS
% - oo_                 [struct]      Dynare's oo_ global structure.
% - errorflag           [logical]     scalar, equal to false iff the simulation did not fail.
%
% REMARKS
% [1] The innovations used for the simulation are saved in oo_.exo_simul, and the resulting paths for the endogenous
%     variables are saved in oo_.endo_simul.
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

debug = false;
errorflag = false;

if ~isempty(innovations)
    oo_.exo_simul(initialconditions.nobs+(1:samplesize),:) = innovations;
end

if ismember(options_.solve_algo, [12,14])
    [funcs, feedback_vars_idxs] = setup_time_recursive_block_simul(M_);
end

function [r, J] = block_wrapper(z, feedback_vars_idx, func, y_dynamic, x, sparse_rowval, sparse_colval, sparse_colptr, T)
    % NB: do as few computations as possible inside this function, since it is
    % called a very large number of times
    y_dynamic(feedback_vars_idx) = z;
    [~, ~, r, J] = feval(func, y_dynamic, x, M_.params, oo_.steady_state, ...
                         sparse_rowval, sparse_colval, sparse_colptr, T);
end

% Simulations (call a Newton-like algorithm for each period).
for it = initialconditions.nobs+(1:samplesize)
    if debug
        dprintf('Period t = %s.', num2str(it-initialconditions.nobs));
    end
    y_ = oo_.endo_simul(:,it-1);
    y = y_;                            % A good guess for the initial conditions is the previous values for the endogenous variables.
    x = oo_.exo_simul(it,:);
    try
        if ismember(options_.solve_algo, [12,14])
            T = NaN(M_.block_structure.dyn_tmp_nbr);
            y_dynamic = [y_; y; NaN(M_.endo_nbr, 1)];
            for blk = 1:length(M_.block_structure.block)
                sparse_rowval = M_.block_structure.block(blk).g1_sparse_rowval;
                sparse_colval = M_.block_structure.block(blk).g1_sparse_colval;
                sparse_colptr = M_.block_structure.block(blk).g1_sparse_colptr;
                if M_.block_structure.block(blk).Simulation_Type ~= 1 % Not an evaluate forward block
                    [z, errorflag, ~, ~, errorcode] = dynare_solve(@block_wrapper, y_dynamic(feedback_vars_idxs{blk}), ...
                                                                   options_.simul.maxit, options_.dynatol.f, ...
                                                                   options_.dynatol.x, options_, ...
                                                                   feedback_vars_idxs{blk}, funcs{blk}, y_dynamic, x, sparse_rowval, sparse_colval, sparse_colptr, T);
                    if errorflag
                        error('Nonlinear solver routine failed with errorcode=%i in block %i and period %i.', errorcode, blk, it)
                    end
                    y_dynamic(feedback_vars_idxs{blk}) = z;
                end
                %% Compute endogenous if the block is of type evaluate or if there are recursive variables in a solve block.
                %% Also update the temporary terms vector.
                [y_dynamic, T] = feval(funcs{blk}, y_dynamic, x, M_.params, ...
                                       oo_.steady_state, sparse_rowval, sparse_colval, ...
                                       sparse_colptr, T);
            end
            oo_.endo_simul(:,it) = y_dynamic(M_.endo_nbr+(1:M_.endo_nbr));
        else
            [oo_.endo_simul(:,it), errorflag, ~, ~, errorcode] = ...
                dynare_solve(@dynamic_backward_model_for_simulation, y, ...
                             options_.simul.maxit, options_.dynatol.f, options_.dynatol.x, ...
                             options_, dynamic_resid, dynamic_g1, y_, x, M_.params, oo_.steady_state, M_.dynamic_g1_sparse_rowval, M_.dynamic_g1_sparse_colval, M_.dynamic_g1_sparse_colptr);
            if errorflag
                error('Nonlinear solver routine failed with errorcode=%i in period %i.', errorcode, it)
            end
        end
    catch Error
        errorflag = true;
        oo_.endo_simul = oo_.endo_simul(:, 1:it-1);
        dprintf('Newton failed on iteration i = %s.', num2str(it-initialconditions.nobs));
        ytm = oo_.endo_simul(:,end);
        xtt = oo_.exo_simul(it,:);
        skipline()
        dprintf('Values of the endogenous variables before the nonlinear solver failure')
        dprintf('----------------------------------------------------------------------')
        skipline()
        dyntable(options_, '', {'VARIABLES','VALUES'}, M_.endo_names(1:M_.orig_endo_nbr), ytm(1:M_.orig_endo_nbr), [], [], 6)
        skipline()
        dprintf('Values of the exogenous variables before the nonlinear solver failure')
        dprintf('---------------------------------------------------------------------')
        skipline()
        dyntable(options_, '', {'VARIABLES','VALUES'}, M_.exo_names, transpose(oo_.exo_simul(it,:)), [], [], 6)
        skipline(2)
        %
        % Get equation tags if any
        %
        if isfield(M_, 'equations_tags')
            etags = cell(M_.orig_endo_nbr, 1);
            for i = 1:M_.orig_endo_nbr
                equations_tags = M_.equations_tags(cellfun(@(x) isequal(x, i), M_.equations_tags(:,1)), :);
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
            etags = split(int2str(1:M_.orig_endo_nbr), '  ');
        end
        %
        % Evaluate and check the residuals
        %
        r = dynamic_backward_model_for_simulation(ytm, dynamic_resid, dynamic_g1, ytm, x, M_.params, oo_.steady_state, M_.dynamic_g1_sparse_rowval, M_.dynamic_g1_sparse_colval, M_.dynamic_g1_sparse_colptr);
        residuals_evaluating_to_nan = isnan(r);
        residuals_evaluating_to_inf = isinf(r);
        residuals_evaluating_to_complex = ~isreal(r);
        if any(residuals_evaluating_to_nan)
            dprintf('Following equations are evaluating to NaN:')
            skipline()
            display_names_of_problematic_equations(M_, eqtags, residuals_evaluating_to_nan);
            skipline()
        end
        if any(residuals_evaluating_to_inf)
            dprintf('Following equations are evaluating to Inf:')
            skipline()
            display_names_of_problematic_equations(M_, eqtags, residuals_evaluating_to_inf);
            skipline()
        end
        if any(residuals_evaluating_to_complex)
            dprintf('Following equations are evaluating to a complex number:')
            skipline()
            display_names_of_problematic_equations(M_, eqtags, residuals_evaluating_to_complex);
            skipline()
        end
        dprintf('Newton failed in period %s with the following error message:', char(initialconditions.lastdate+(it-initialconditions.nobs)));
        skipline()
        dprintf('\t %s', Error.message);
        skipline()
        break
        % TODO Implement same checks with the jacobian matrix.
        % TODO Modify other solvers to return an exitflag.
    end
end

ysim = oo_.endo_simul(1:M_.orig_endo_nbr,:);
xsim = oo_.exo_simul;

end

function display_names_of_problematic_equations(M_, eqtags, TruthTable)
for i=1:M_.orig_endo_nbr
    if TruthTable(i)
        dprintf(' - %s', eqtags{i})
    end
end
for i=M_.orig_endo_nbr+1:M_.endo_nbr
    if TruthTable(i)
        dprintf(' - Auxiliary equation for %s', M_.endo_names{i})
    end
end
end
