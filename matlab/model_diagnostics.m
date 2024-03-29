function model_diagnostics(M_,options_,oo_)
% function model_diagnostics(M_,options_,oo_)
%   computes various diagnostics on the model
% INPUTS
%   M_         [matlab structure] Definition of the model.
%   options_   [matlab structure] Global options.
%   oo_        [matlab structure] Results
%
% OUTPUTS
%   none
%
% ALGORITHM
%   ...
%
% SPECIAL REQUIREMENTS
%   none.
%

% Copyright © 1996-2024 Dynare Team
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

endo_names = M_.endo_names;
lead_lag_incidence = M_.lead_lag_incidence;
maximum_endo_lag = M_.maximum_endo_lag;

if options_.ramsey_policy
    %test whether specification matches
    inst_nbr = size(options_.instruments,1);
    if inst_nbr~=0
        implied_inst_nbr = M_.ramsey_orig_endo_nbr - M_.ramsey_orig_eq_nbr;
        if inst_nbr>implied_inst_nbr
            warning('You have specified more steady state instruments than there are omitted equations. While there are use cases for this setup, it is rather unusual. Check whether this is desired.')
        elseif inst_nbr<implied_inst_nbr
            warning('You have specified fewer steady state instruments than there are omitted equations. While there are use cases for this setup, it is rather unusual. Check whether this is desired.')
        end
    else
        if options_.steadystate_flag
            warning('You have specified a steady state file, but not provided steady state instruments. In this case, you typically need to make sure to provide all steady state values, including the ones for the planner''s instrument(s).')
        end
    end
end

problem_dummy=0;

%naming conflict in steady state file
if options_.steadystate_flag == 1
    if strmatch('ys',M_.endo_names,'exact') 
        disp('MODEL_DIAGNOSTICS: using the name ys for an endogenous variable will typically conflict with the internal naming in user-defined steady state files.')
        problem_dummy=1;
    end
    if strmatch('ys',M_.param_names,'exact')
        disp('MODEL_DIAGNOSTICS: using the name ys for a parameter will typically conflict with the internal naming in user-defined steady state files.')
        problem_dummy=1;
    end
    if strmatch('M_',M_.endo_names,'exact') 
        disp('MODEL_DIAGNOSTICS: using the name M_ for an endogenous variable will typically conflict with the internal naming in user-defined steady state files.')
        problem_dummy=1;
    end
    if strmatch('M_',M_.param_names,'exact')
        disp('MODEL_DIAGNOSTICS: using the name M_ for a parameter will typically conflict with the internal naming in user-defined steady state files.')
        problem_dummy=1;
    end
end

%
% missing variables at the current period
%
k = find(lead_lag_incidence(maximum_endo_lag+1,:)==0);
if ~isempty(k)
    problem_dummy=1;
    disp(['MODEL_DIAGNOSTICS: The following endogenous variables aren''t present at ' ...
          'the current period in the model:'])
    for i=1:length(k)
        disp(endo_names{k(i)})
    end
end

%
% check steady state
%
info = 0;

if M_.exo_nbr == 0
    oo_.exo_steady_state = [] ;
end


info=test_for_deep_parameters_calibration(M_);
if info
    problem_dummy=1;
end

% check if ys is steady state
options_.debug=true; %locally set debug option to true
if options_.logged_steady_state %if steady state was previously logged, undo this
    oo_.dr.ys=exp(oo_.dr.ys);
    oo_.steady_state=exp(oo_.steady_state);
    options_.logged_steady_state=0;
end
[dr.ys,M_.params,check1]=evaluate_steady_state(oo_.steady_state,[oo_.exo_steady_state; oo_.exo_det_steady_state],M_,options_,~options_.steadystate.nocheck);

if isfield(M_,'occbin')
    if any(oo_.exo_steady_state)
        disp('MODEL_DIAGNOSTICS: OccBin was detected in conjunction with a non-zero steady state of the exogenous variables. That will usually create issues.')
        problem_dummy=1;
    end
end
% testing for problem
if check1(1)
    problem_dummy=1;
    disp('MODEL_DIAGNOSTICS: The steady state cannot be computed')
    if any(isnan(dr.ys))
        disp('MODEL_DIAGNOSTICS: Steady state contains NaNs')
    end
    if any(isinf(dr.ys))
        disp('MODEL_DIAGNOSTICS: Steady state contains Inf')
    end
    return
end

if ~isreal(dr.ys)
    problem_dummy=1;
    disp(['MODEL_DIAGNOSTICS: Steady state contains complex ' ...
          'numbers'])
    return
end

%
% singular Jacobian of static model
%
singularity_problem = 0;
if ~options_.block
    nb = 1;
else
    nb = length(M_.block_structure_stat.block);
end

exo = [oo_.exo_steady_state; oo_.exo_det_steady_state];
for b=1:nb
    if options_.bytecode
        if nb == 1
            [~, jacob] = bytecode(M_, options_, dr.ys, exo, M_.params, dr.ys, 1, exo, ...
                                    'evaluate', 'static');
        else
            [~, jacob] = bytecode(M_, options_, dr.ys, exo, M_.params, dr.ys, 1, exo, ...
                                    'evaluate', 'static', 'block_decomposed', ['block=' ...
                                int2str(b)]);
        end
    else
        if options_.block
            T = NaN(M_.block_structure_stat.tmp_nbr, 1);
            fh_static = str2func(sprintf('%s.sparse.block.static_%d', M_.fname, b));
            [~, ~,~, jacob] = fh_static(dr.ys, exo, M_.params, M_.block_structure_stat.block(b).g1_sparse_rowval, ...
                M_.block_structure_stat.block(b).g1_sparse_colval, ...
                M_.block_structure_stat.block(b).g1_sparse_colptr, T);
            n_vars_jacob=size(jacob,2);
        else
            [~, T_order, T] = feval([M_.fname '.sparse.static_resid'], dr.ys, exo, M_.params);
            jacob = feval([M_.fname '.sparse.static_g1'], dr.ys, exo, M_.params, M_.static_g1_sparse_rowval, M_.static_g1_sparse_colval, M_.static_g1_sparse_colptr, T_order, T);
            n_vars_jacob=M_.endo_nbr;
        end
        jacob=full(jacob);
    end
    if any(any(isinf(jacob) | isnan(jacob)))
        problem_dummy=1;
        [infrow,infcol]=find(isinf(jacob) | isnan(jacob));
        fprintf('\nMODEL_DIAGNOSTICS: The Jacobian of the static model contains Inf or NaN. The problem arises from: \n\n')
        display_problematic_vars_Jacobian(infrow,infcol,M_,dr.ys,'static','MODEL_DIAGNOSTICS: ')
    end
    if any(any(~isreal(jacob)))
        problem_dummy=1;
        [imagrow,imagcol]=find(abs(imag(jacob))>1e-15);
        fprintf('\nMODEL_DIAGNOSTICS: The Jacobian of the static model contains imaginary parts. The problem arises from: \n\n')
        display_problematic_vars_Jacobian(imagrow,imagcol,M_,dr.ys,'static','MODEL_DIAGNOSTICS: ')
    end
    try
        if (~isoctave && matlab_ver_less_than('9.12')) || isempty(options_.jacobian_tolerance)
            rank_jacob = rank(jacob); %can sometimes fail
        else
            rank_jacob = rank(jacob,options_.jacobian_tolerance); %can sometimes fail
        end
    catch
        rank_jacob=size(jacob,1);
    end
    if rank_jacob < size(jacob,1)
        problem_dummy=1;
        singularity_problem = 1;
        disp(['MODEL_DIAGNOSTICS:  The Jacobian of the static model is ' ...
              'singular'])
        disp(['MODEL_DIAGNOSTICS:  there is ' num2str(n_vars_jacob-rank_jacob) ...
              ' collinear relationships between the variables and the equations'])
        if (~isoctave && matlab_ver_less_than('9.12')) || isempty(options_.jacobian_tolerance)
            ncol = null(jacob);
        else
            ncol = null(jacob,options_.jacobian_tolerance); %can sometimes fail
        end
        n_rel = size(ncol,2);
        for i = 1:n_rel
            if n_rel  > 1
                disp(['Relation ' int2str(i)])
            end
            disp('Collinear variables:')
            for j=1:10
                k = find(abs(ncol(:,i)) > 10^-j);
                if max(abs(jacob(:,k)*ncol(k,i))) < 1e-6
                    break
                end
            end
            if options_.block && ~options_.bytecode
                fprintf('%s\n',endo_names{M_.block_structure_stat.block(b).variable(k)})
            else
                fprintf('%s\n',endo_names{k})
            end
        end
        if (~isoctave && matlab_ver_less_than('9.12')) || isempty(options_.jacobian_tolerance)
            neq = null(jacob'); %can sometimes fail
        else
            neq = null(jacob',options_.jacobian_tolerance); %can sometimes fail
        end
        n_rel = size(neq,2);
        for i = 1:n_rel
            if n_rel  > 1
                disp(['Relation ' int2str(i)])
            end
            disp('Collinear equations')
            for j=1:10
                k = find(abs(neq(:,i)) > 10^-j);
                if max(abs(jacob(k,:)'*neq(k,i))) < 1e-6
                    break
                end
            end
            if options_.block && ~options_.bytecode
                disp(M_.block_structure_stat.block(b).equation(k))
            else
                disp(k')
            end
        end
    end
end

if singularity_problem
    try
        options_check=options_;
        options_check.noprint=1;
        [eigenvalues_] = check(M_, options_check, oo_);
        if any(abs(abs(eigenvalues_)-1)<1e-6)
            fprintf('MODEL_DIAGNOSTICS:  The singularity seems to be (partly) caused by the presence of a unit root\n')
            fprintf('MODEL_DIAGNOSTICS:  as the absolute value of one eigenvalue is in the range of +-1e-6 to 1.\n')
            fprintf('MODEL_DIAGNOSTICS:  If the model is actually supposed to feature unit root behavior, such a warning is expected,\n')
            fprintf('MODEL_DIAGNOSTICS:  but you should nevertheless check whether there is an additional singularity problem.\n')
        end
    catch
    end
    fprintf('MODEL_DIAGNOSTICS:  The presence of a singularity problem typically indicates that there is one\n')
    fprintf('MODEL_DIAGNOSTICS:  redundant equation entered in the model block, while another non-redundant equation\n')
    fprintf('MODEL_DIAGNOSTICS:  is missing. The problem often derives from Walras Law.\n')
end

%%check dynamic Jacobian
dyn_endo_ss = repmat(dr.ys, 3, 1);

if options_.order == 1
    if (options_.bytecode)
        [~, loc_dr] = bytecode('dynamic','evaluate', M_, options_, z, exo_simul, ...
                               M_.params, dr.ys, 1);
        % TODO: simplify the following once bytecode MEX has been updated to sparse format
        g1 = zeros(M_.endo_nbr, 3*M_.endo_nbr+M_.exo_nbr+M_.exo_det_nbr);
        if M_.maximum_endo_lag > 0
            g1(:, find(M_.lead_lag_incidence(M_.maximum_endo_lag, :))) = loc_dr.g1(:, 1:M_.nspred);
        end
        [~,icurr] = find(M_.lead_lag_incidence(M_.maximum_endo_lag+1, :));
        g1(:, M_.endo_nbr + icurr) = loc_dr.g1(:, M_.nspred+(1:length(icurr)));
        if M_.maximum_endo_lead > 0
            g1(:, 2*M_.endo_nbr + find(M_.lead_lag_incidence(M_.maximum_endo_lag+2, :))) = loc_dr.g1(:, M_.nspred+M_.endo_nbr+(1:M_.nsfwrd));
        end
        g1(:, 3*M_.endo_nbr+(1:M_.exo_nbr)) = loc_dr.g1_x;
        g1(:, 3*M_.endo_nbr+M_.exo_nbr+(1:M_.exo_det_nbr)) = loc_dr.g1_xd;
        g1 = sparse(g1);
    else
        g1 = feval([M_.fname '.sparse.dynamic_g1'], dyn_endo_ss, exo, M_.params, dr.ys, ...
                   M_.dynamic_g1_sparse_rowval, M_.dynamic_g1_sparse_colval, ...
                   M_.dynamic_g1_sparse_colptr);
    end
elseif options_.order >= 2
    [g1, T_order, T] = feval([M_.fname '.sparse.dynamic_g1'], dyn_endo_ss, exo, M_.params, ...
                             dr.ys, M_.dynamic_g1_sparse_rowval, M_.dynamic_g1_sparse_colval, ...
                             M_.dynamic_g1_sparse_colptr);
    g2_v = feval([M_.fname '.sparse.dynamic_g2'], dyn_endo_ss, exo, M_.params, dr.ys, T_order, T);
end

if any(any(isinf(g1) | isnan(g1)))
    problem_dummy=1;
    [infrow,infcol]=find(isinf(g1) | isnan(g1));
    fprintf('\nMODEL_DIAGNOSTICS: The Jacobian of the dynamic model contains Inf or NaN. The problem arises from: \n\n')
    display_problematic_vars_Jacobian(infrow,infcol,M_,dr.ys,'dynamic','MODEL_DIAGNOSTICS: ')
end
if any(any(~isreal(g1)))
    [imagrow,imagcol]=find(abs(imag(g1))>1e-15);
    if ~isempty(imagrow)
        problem_dummy=1;
        fprintf('\nMODEL_DIAGNOSTICS: The Jacobian of the dynamic model contains imaginary parts. The problem arises from: \n\n')
        display_problematic_vars_Jacobian(imagrow,imagcol,M_,dr.ys,'dynamic','MODEL_DIAGNOSTICS: ')
    end
end
if exist('g2_v','var')
    if any(any(isinf(g2_v) | isnan(g2_v)))
        problem_dummy=1;
        fprintf('\nMODEL_DIAGNOSTICS: The Hessian of the dynamic model contains Inf or NaN.\n')
    end
end

if problem_dummy==0
    fprintf('MODEL_DIAGNOSTICS:  No obvious problems with this mod-file were detected.\n')
end
