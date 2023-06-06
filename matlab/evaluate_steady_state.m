function [ys,params,info] = evaluate_steady_state(ys_init,M,options,oo,steadystate_check_flag)
% function [ys,params,info] = evaluate_steady_state(ys_init,M,options,oo,steadystate_check_flag)
% Computes the steady state
%
% INPUTS
%   ys_init                   vector           initial values used to compute the steady
%                                                 state
%   M                         struct           model structure
%   options                   struct           options
%   oo                        struct           output results
%   steadystate_check_flag    boolean          if true, check that the
%                                              steadystate verifies the
%                                              static model
%
% OUTPUTS
%   ys                        vector           steady state (in declaration order)
%   params                    vector           model parameters possibly
%                                              modified by user steadystate
%                                              function
%   info                      2x1 vector       error codes
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2001-2023 Dynare Team
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

if options.solve_algo < 0 || options.solve_algo > 14
    error('STEADY: solve_algo must be between 0 and 14')
end

if ~options.bytecode && ~options.block && options.solve_algo > 4 && ...
        options.solve_algo < 9
    error('STEADY: you can''t use solve_algo = {5,6,7,8} without block nor bytecode options')
end

if ~options.bytecode && options.block && options.solve_algo == 5
    error('STEADY: you can''t use solve_algo = 5 without bytecode option')
end

if isoctave && options.solve_algo == 11
    error(['STEADY: you can''t use solve_algo = %u under Octave'],options.solve_algo)
end

info = 0;
check = 0;

steadystate_flag = options.steadystate_flag;
params = M.params;
exo_ss = [oo.exo_steady_state; oo.exo_det_steady_state];

if length(M.aux_vars) > 0 && ~steadystate_flag && M.set_auxiliary_variables
    h_set_auxiliary_variables = str2func([M.fname '.set_auxiliary_variables']);
    ys_init = h_set_auxiliary_variables(ys_init,exo_ss,params);
end

if options.ramsey_policy
    if ~isfinite(M.params(strmatch('optimal_policy_discount_factor',M.param_names,'exact')))
        fprintf('\nevaluate_steady_state: the planner_discount is NaN/Inf. That will cause problems.\n')
    end
    if steadystate_flag
        % explicit steady state file
        [ys,params,info] = evaluate_steady_state_file(ys_init,exo_ss,M, ...
                                                      options,steadystate_check_flag);
        %test whether it solves model conditional on the instruments
        if ~options.debug
            resids = evaluate_static_model(ys,exo_ss,params,M,options);
        else
            [resids, ~ , jacob]= evaluate_static_model(ys,exo_ss,params,M,options);
        end
        nan_indices=find(isnan(resids(M.ramsey_orig_endo_nbr+(1:M.ramsey_orig_eq_nbr))));

        if ~isempty(nan_indices)
            if options.debug
                fprintf('\nevaluate_steady_state: The steady state file computation for the Ramsey problem resulted in NaNs.\n')
                fprintf('evaluate_steady_state: The steady state was computed conditional on the following initial instrument values: \n')
                for ii = 1:size(options.instruments,1)
                    fprintf('\t %s \t %f \n',options.instruments{ii},ys_init(strmatch(options.instruments{ii},M.endo_names,'exact')))
                end
                fprintf('evaluate_steady_state: The problem occured in the following equations: \n')
                fprintf('\t Equation(s): ')
                for ii=1:length(nan_indices)
                    fprintf('%d, ',nan_indices(ii));
                end
                skipline()
                fprintf('evaluate_steady_state: If those initial values are not admissable, change them using an initval-block.\n')
                skipline(2)
            end
            info(1) = 84;
            info(2) = resids'*resids;
            return
        end

        if any(imag(ys(M.ramsey_orig_endo_nbr+(1:M.ramsey_orig_eq_nbr))))
            if options.debug
                fprintf('\nevaluate_steady_state: The steady state file computation for the Ramsey problem resulted in complex numbers.\n')
                fprintf('evaluate_steady_state: The steady state was computed conditional on the following initial instrument values: \n')
                for ii = 1:size(options.instruments,1)
                    fprintf('\t %s \t %f \n',options.instruments{ii},ys_init(strmatch(options.instruments{ii},M.endo_names,'exact')))
                end
                fprintf('evaluate_steady_state: If those initial values are not admissable, change them using an initval-block.\n')
                skipline(2)
            end
            info(1) = 86;
            info(2) = resids'*resids;
            return
        end

        if max(abs(resids(M.ramsey_orig_endo_nbr+(1:M.ramsey_orig_eq_nbr)))) > options.solve_tolf %does it solve for all variables except for the Lagrange multipliers
            if options.debug
                fprintf('\nevaluate_steady_state: The steady state file does not solve the steady state for the Ramsey problem.\n')
                fprintf('evaluate_steady_state: Conditional on the following instrument values: \n')
                for ii = 1:size(options.instruments,1)
                    fprintf('\t %s \t %f \n',options.instruments{ii},ys_init(strmatch(options.instruments{ii},M.endo_names,'exact')))
                end
                fprintf('evaluate_steady_state: the following equations have non-zero residuals: \n')
                for ii=M.ramsey_orig_endo_nbr+1:M.endo_nbr
                    if abs(resids(ii)) > options.solve_tolf
                        fprintf('\t Equation number %d: %f\n',ii-M.ramsey_orig_endo_nbr, resids(ii))
                    end
                end
                skipline(2)
            end
            info(1) = 85;
            info(2) = resids'*resids;
            return
        end
    end
    if options.debug
        if steadystate_flag
            infrow=find(isinf(ys_init(1:M.orig_endo_nbr)));
        else
            infrow=find(isinf(ys_init));
        end
        if ~isempty(infrow)
            fprintf('\nevaluate_steady_state: The initial values for the steady state of the following variables are Inf:\n');
            for iter=1:length(infrow)
                fprintf('%s\n',M.endo_names{infrow(iter)});
            end
        end
        if steadystate_flag
            nanrow=find(isnan(ys_init(1:M.orig_endo_nbr)));
        else
            nanrow=find(isnan(ys_init));
        end
        if ~isempty(nanrow)
            fprintf('\nevaluate_steady_state: The initial values for the steady state of the following variables are NaN:\n');
            for iter=1:length(nanrow)
                fprintf('%s\n',M.endo_names{nanrow(iter)});
            end
        end
        if steadystate_flag
            nan_indices_mult=find(isnan(resids(1:M.ramsey_orig_endo_nbr)));
            if any(nan_indices_mult)
                fprintf('evaluate_steady_state: The steady state results NaN for auxiliary equation %u.\n',nan_indices_mult);
                fprintf('evaluate_steady_state: This is often a sign of problems.\n');
            end
            [infrow,infcol]=find(isinf(jacob));
            
            if ~isempty(infrow)
                fprintf('\nevaluate_steady_state: The Jacobian of the dynamic model contains Inf. The problem is associated with:\n\n')
                display_problematic_vars_Jacobian(infrow,infcol,M,ys,'static','evaluate_steady_state: ')
            end
            
            if ~isreal(jacob)
                [imagrow,imagcol]=find(abs(imag(jacob))>1e-15);
                fprintf('\nevaluate_steady_state: The Jacobian of the dynamic model contains imaginary parts. The problem arises from: \n\n')
                display_problematic_vars_Jacobian(imagrow,imagcol,M,ys,'static','evaluate_steady_state: ')
            end
            
            [nanrow,nancol]=find(isnan(jacob));
            if ~isempty(nanrow)
                fprintf('\nevaluate_steady_state: The Jacobian of the dynamic model contains NaN. The problem is associated with:\n\n')
                display_problematic_vars_Jacobian(nanrow,nancol,M,ys,'static','evaluate_steady_state: ')
            end
            
        end

    end
    %either if no steady state file or steady state file without problems
    [ys,params,info] = dyn_ramsey_static(ys_init,M,options,oo);
    if info
        return
    end
    %check whether steady state really solves the model
    resids = evaluate_static_model(ys,exo_ss,params,M,options);

    nan_indices_multiplier=find(isnan(resids(1:M.ramsey_orig_endo_nbr)));
    nan_indices=find(isnan(resids(M.ramsey_orig_endo_nbr+1:end)));

    if ~isempty(nan_indices)
        if options.debug
            fprintf('\nevaluate_steady_state: The steady state computation for the Ramsey problem resulted in NaNs.\n')
            fprintf('evaluate_steady_state: The steady state computation resulted in the following instrument values: \n')
            for i = 1:size(options.instruments,1)
                fprintf('\t %s \t %f \n',options.instruments{i},ys(strmatch(options.instruments{i},M.endo_names,'exact')))
            end
            fprintf('evaluate_steady_state: The problem occured in the following equations: \n')
            fprintf('\t Equation(s): ')
            for ii=1:length(nan_indices)
                fprintf('%d, ',nan_indices(ii));
            end
            skipline()
        end
        info(1) = 82;
        return
    end

    if ~isempty(nan_indices_multiplier)
        if options.debug
            fprintf('\nevaluate_steady_state: The steady state computation for the Ramsey problem resulted in NaNs in the auxiliary equations.\n')
            fprintf('evaluate_steady_state: The steady state computation resulted in the following instrument values: \n')
            for i = 1:size(options.instruments,1)
                fprintf('\t %s \t %f \n',options.instruments{i},ys(strmatch(options.instruments{i},M.endo_names,'exact')))
            end
            fprintf('evaluate_steady_state: The problem occured in the following equations: \n')
            fprintf('\t Auxiliary equation(s): ')
            for ii=1:length(nan_indices_multiplier)
                fprintf('%d, ',nan_indices_multiplier(ii));
            end
            skipline()
        end
        info(1) = 83;
        return
    end

    if max(abs(resids)) > options.solve_tolf %does it solve for all variables including the auxiliary ones
        if options.debug
            fprintf('\nevaluate_steady_state: The steady state for the Ramsey problem could not be computed.\n')
            fprintf('evaluate_steady_state: The steady state computation stopped with the following instrument values:: \n')
            for i = 1:size(options.instruments,1)
                fprintf('\t %s \t %f \n',options.instruments{i},ys(strmatch(options.instruments{i},M.endo_names,'exact')))
            end
            fprintf('evaluate_steady_state: The following equations have non-zero residuals: \n')
            for ii=1:M.ramsey_orig_endo_nbr
                if abs(resids(ii)) > options.solve_tolf/100
                    fprintf('\t Auxiliary Ramsey equation number %d: %f\n',ii, resids(ii))
                end
            end
            for ii=M.ramsey_orig_endo_nbr+1:M.endo_nbr
                if abs(resids(ii)) > options.solve_tolf/100
                    fprintf('\t Equation number %d: %f\n',ii-M.ramsey_orig_endo_nbr, resids(ii))
                end
            end
            skipline(2)
        end
        info(1) = 81;
        info(2) = resids'*resids;
        return
    end
elseif steadystate_flag
    % explicit steady state file
    [ys,params,info] = evaluate_steady_state_file(ys_init,exo_ss,M, options,steadystate_check_flag);
    if size(ys,2)>size(ys,1)
        error('STEADY: steady_state-file must return a column vector, not a row vector.')
    end
    if info(1)
        return
    end
elseif ~options.bytecode && ~options.block
    static_resid = str2func(sprintf('%s.sparse.static_resid', M.fname));
    static_g1 = str2func(sprintf('%s.sparse.static_g1', M.fname));
    if ~options.linear
        % non linear model
        if  ismember(options.solve_algo,[10,11])
            [lb,ub,eq_index] = get_complementarity_conditions(M,options.ramsey_policy);
            if options.solve_algo == 10
                options.lmmcp.lb = lb;
                options.lmmcp.ub = ub;
            elseif options.solve_algo == 11
                options.mcppath.lb = lb;
                options.mcppath.ub = ub;
            end
            [ys,check,fvec] = dynare_solve(@static_mcp_problem,...
                ys_init,...
                options.steady.maxit, options.solve_tolf, options.solve_tolx, ...
                options, exo_ss, params,...
                M.endo_nbr, static_resid, static_g1, ...
                M.static_g1_sparse_rowval, M.static_g1_sparse_colval, M.static_g1_sparse_colptr, eq_index);
        else
            [ys, check] = dynare_solve(@static_problem, ys_init, ...
                options.steady.maxit, options.solve_tolf, options.solve_tolx, ...
                options, exo_ss, params, M.endo_nbr, static_resid, static_g1, ...
                M.static_g1_sparse_rowval, M.static_g1_sparse_colval, M.static_g1_sparse_colptr);
        end
        if check && options.debug
            [ys, check, fvec, fjac, errorcode] = dynare_solve(@static_problem, ys_init, ...
                                                              options.steady.maxit, options.solve_tolf, options.solve_tolx, ...
                                                              options, exo_ss, params, M.endo_nbr, static_resid, static_g1, M.static_g1_sparse_rowval, M.static_g1_sparse_colval, M.static_g1_sparse_colptr);
            dprintf('Nonlinear solver routine returned errorcode=%i.', errorcode)
            skipline()
            [infrow,infcol]=find(isinf(fjac) | isnan(fjac));
            if ~isempty(infrow)
                fprintf('\nSTEADY:  The Jacobian at the initial values contains Inf or NaN. The problem arises from: \n')
                display_problematic_vars_Jacobian(infrow,infcol,M,ys_init,'static','STEADY: ')
            end
            problematic_equation = find(~isfinite(fvec));
            if ~isempty(problematic_equation)
                fprintf('\nSTEADY:  numerical initial values or parameters incompatible with the following equations\n')
                disp(problematic_equation')
                fprintf('Please check for example\n')
                fprintf('   i) if all parameters occurring in these equations are defined\n')
                fprintf('  ii) that no division by an endogenous variable initialized to 0 occurs\n')
            end
        end
    else
        % linear model
        [fvec, T_order, T] = static_resid(ys_init, exo_ss, params);
        jacob = static_g1(ys_init, exo_ss, params, M.static_g1_sparse_rowval, M.static_g1_sparse_colval, M.static_g1_sparse_colptr, T_order, T);

        ii = find(~isfinite(fvec));
        if ~isempty(ii)
            ys=fvec;
            check=1;
            disp(['STEADY:  numerical initial values or parameters incompatible with the following' ...
                  ' equations'])
            disp(ii')
            disp('Check whether your model is truly linear. Put "resid(1);" before "steady;" to see the problematic equations.')
        elseif isempty(ii) && max(abs(fvec)) > 1e-12
            ys = ys_init-jacob\fvec;
            resid = evaluate_static_model(ys,exo_ss,params,M,options);
            if max(abs(resid)) > 1e-6
                check=1;
                fprintf('STEADY: No steady state for your model could be found\n')
                fprintf('STEADY: Check whether your model is truly linear. Put "resid(1);" before "steady;" to see the problematic equations.\n')
            end
        else
            ys = ys_init;
        end
        if options.debug
            if any(any(isinf(jacob) | isnan(jacob)))
                [infrow,infcol]=find(isinf(jacob) | isnan(jacob));
                fprintf('\nSTEADY:  The Jacobian contains Inf or NaN. The problem arises from: \n\n')
                for ii=1:length(infrow)
                    fprintf('STEADY:  Derivative of Equation %d with respect to Variable %s  (initial value of %s: %g) \n',infrow(ii),M.endo_names{infcol(ii),:},M.endo_names{infcol(ii),:},ys_init(infcol(ii)))
                end
                fprintf('Check whether your model is truly linear. Put "resid(1);" before "steady;" to see the problematic equations.\n')
            end
        end
    end
elseif ~options.bytecode && options.block
    ys = ys_init;
    T = NaN(M.block_structure_stat.tmp_nbr, 1);
    for b = 1:length(M.block_structure_stat.block)
        fh_static = str2func(sprintf('%s.sparse.block.static_%d', M.fname, b));
        if M.block_structure_stat.block(b).Simulation_Type ~= 1 && ...
                M.block_structure_stat.block(b).Simulation_Type ~= 2
            mfs_idx = M.block_structure_stat.block(b).variable(end-M.block_structure_stat.block(b).mfs+1:end);
            if options.solve_algo <= 4 || options.solve_algo >= 9
                [ys(mfs_idx), errorflag] = dynare_solve(@block_mfs_steadystate, ys(mfs_idx), ...
                                                        options.simul.maxit, options.solve_tolf, options.solve_tolx, ...
                                                        options, fh_static, b, ys, exo_ss, params, T, M);
                if errorflag
                    check = 1;
                    break
                end
            else
                nze = length(M.block_structure_stat.block(b).g1_sparse_rowval);
                [ys, T, success] = solve_one_boundary(fh_static, ys, exo_ss, ...
                                                      params, [], T, mfs_idx, nze, 1, false, b, 0, options.simul.maxit, ...
                                                      options.solve_tolf, ...
                                                      0, options.solve_algo, true, false, false, M, options);
                if ~success
                    check = 1;
                    break
                end
            end
        end
        % Compute endogenous if the block is of type evaluate forward/backward or if there are recursive variables in a solve block.
        % Also update the temporary terms vector (needed for the dynare_solve case)
        [ys, T] = fh_static(ys, exo_ss, params, M.block_structure_stat.block(b).g1_sparse_rowval, ...
                            M.block_structure_stat.block(b).g1_sparse_colval, ...
                            M.block_structure_stat.block(b).g1_sparse_colptr, T);
    end
elseif options.bytecode
    if options.solve_algo >= 5 && options.solve_algo <= 8
        try
            if options.block
                ys = bytecode('static', 'block_decomposed', ys_init, exo_ss, params);
            else
                ys = bytecode('static', ys_init, exo_ss, params);
            end
        catch ME
            disp(ME.message);
            check = 1;
        end
    elseif options.block
        ys = ys_init;
        T = NaN(M.block_structure_stat.tmp_nbr, 1);
        for b = 1:length(M.block_structure_stat.block)
            if M.block_structure_stat.block(b).Simulation_Type ~= 1 && ...
                    M.block_structure_stat.block(b).Simulation_Type ~= 2
                mfs_idx = M.block_structure_stat.block(b).variable(end-M.block_structure_stat.block(b).mfs+1:end);
                [ys(mfs_idx), errorflag] = dynare_solve(@block_bytecode_mfs_steadystate, ...
                                                        ys(mfs_idx), options.simul.maxit, ...
                                                        options.solve_tolf, options.solve_tolx, ...
                                                        options, b, ys, exo_ss, params, T, M);
                if errorflag
                    check = 1;
                    break
                end
            end
            % Compute endogenous if the block is of type evaluate forward/backward or if there are recursive variables in a solve block.
            % Also update the temporary terms vector (needed for the dynare_solve case)
            try
                [~, ~, ys, T] = bytecode(ys, exo_ss, params, ys, 1, ys, T, 'evaluate', 'static', ...
                                         'block_decomposed', ['block = ' int2str(b)]);
            catch ME
                disp(ME.message);
                check = 1;
                break
            end
        end
    else
        [ys, check] = dynare_solve(@bytecode_steadystate, ys_init, ...
                                   options.simul.maxit, options.solve_tolf, options.solve_tolx, ...
                                   options, exo_ss, params);
    end
end

if check
    info(1)= 20;
    %make sure ys contains auxiliary variables in case of problem with dynare_solve
    if length(M.aux_vars) > 0 && ~steadystate_flag
        if M.set_auxiliary_variables
            ys = h_set_auxiliary_variables(ys,exo_ss,params);
        end
    end
    resid = evaluate_static_model(ys,exo_ss,params,M,options);
    info(2) = resid'*resid ;
    if isnan(info(2))
        info(1)=22;
    end
    return
end

% If some equations are tagged [static] or [dynamic], verify consistency
if M.static_and_dynamic_models_differ
    % Evaluate residual of *dynamic* model using the steady state
    % computed on the *static* one
    if options.bytecode
        z = repmat(ys,1,M.maximum_lead + M.maximum_lag + 1);
        zx = repmat([exo_ss'], M.maximum_lead + M.maximum_lag + 1, 1);
        [r, ~]= bytecode('dynamic','evaluate', z, zx, params, ys, 1);
    else
        r = feval([M.fname '.sparse.dynamic_resid'], repmat(ys, 3, 1), exo_ss, params, ys);
    end
    % Fail if residual greater than tolerance
    if max(abs(r)) > options.solve_tolf
        info(1) = 25;
        return
    end
end

if ~isreal(ys)
    info(1) = 21;
    info(2) = sum(imag(ys).^2);
    ys = real(ys);
    return
end

if ~isempty(find(isnan(ys)))
    info(1) = 22;
    info(2) = NaN;
    return
end

function [resids,jac] = static_problem(y, x, params, nvar, fh_static_resid, fh_static_g1, sparse_rowval, sparse_colval, sparse_colptr)
[r, T_order, T] = fh_static_resid(y, x, params);
j = fh_static_g1(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T_order, T);
resids = r(1:nvar);
jac = j(1:nvar,1:nvar);

function [resids,jac] = static_mcp_problem(y, x, params, nvar, fh_static_resid, fh_static_g1, sparse_rowval, sparse_colval, sparse_colptr, eq_index)
[r, T_order, T] = fh_static_resid(y, x, params);
j = fh_static_g1(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T_order, T);
resids = r(eq_index);
jac = j(eq_index,1:nvar);

function [r, g1] = block_mfs_steadystate(y, fh_static, b, y_all, exo, params, T, M)
% Wrapper around the static files, for block without bytecode
mfs_idx = M.block_structure_stat.block(b).variable(end-M.block_structure_stat.block(b).mfs+1:end);
y_all(mfs_idx) = y;
[~,~,r,g1] = fh_static(y_all, exo, params, M.block_structure_stat.block(b).g1_sparse_rowval, ...
                       M.block_structure_stat.block(b).g1_sparse_colval, ...
                       M.block_structure_stat.block(b).g1_sparse_colptr, T);

function [r, g1] = bytecode_steadystate(y, exo, params)
% Wrapper around the static file, for bytecode (without block)
[r, g1] = bytecode(y, exo, params, y, 1, exo, 'evaluate', 'static');

function [r, g1] = block_bytecode_mfs_steadystate(y, b, y_all, exo, params, T, M)
% Wrapper around the static files, for block without bytecode
mfs_idx = M.block_structure_stat.block(b).variable(end-M.block_structure_stat.block(b).mfs+1:end);
y_all(mfs_idx) = y;
[r, g1] = bytecode(y_all, exo, params, y_all, 1, y_all, T, 'evaluate', 'static', 'block_decomposed', ['block = ' int2str(b) ]);
g1 = g1(:,end-M.block_structure_stat.block(b).mfs+1:end); % Make Jacobian square if mfs>0
