function model_diagnostics(M,options,oo)
% function model_diagnostics(M,options,oo)
%   computes various diagnostics on the model
% INPUTS
%   M         [matlab structure] Definition of the model.
%   options   [matlab structure] Global options.
%   oo        [matlab structure] Results
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

% Copyright © 1996-2020 Dynare Team
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

endo_nbr = M.endo_nbr;
endo_names = M.endo_names;
lead_lag_incidence = M.lead_lag_incidence;
maximum_endo_lag = M.maximum_endo_lag;

if options.ramsey_policy
    %test whether specification matches
    inst_nbr = size(options.instruments,1);
    if inst_nbr~=0
        orig_endo_aux_nbr = M.orig_endo_nbr + min(find([M.aux_vars.type] == 6)) - 1;
        implied_inst_nbr = orig_endo_aux_nbr - M.orig_eq_nbr;
        if inst_nbr>implied_inst_nbr
            error('You have specified more instruments than there are omitted equations')
        elseif inst_nbr<implied_inst_nbr
            error('You have specified fewer instruments than there are omitted equations')
        end
    else
        if options.steadystate_flag
            error('You have specified a steady state file, but not provided an instrument. Either delete the steady state file or provide an instrument')
        end
    end
end

problem_dummy=0;

%naming conflict in steady state file
if options.steadystate_flag == 1
    if strmatch('ys',M.endo_names,'exact') 
        disp(['MODEL_DIAGNOSTICS: using the name ys for an endogenous variable will typically conflict with the internal naming in user-defined steady state files.'])
        problem_dummy=1;
    end
    if strmatch('ys',M.param_names,'exact')
        disp(['MODEL_DIAGNOSTICS: using the name ys for a parameter will typically conflict with the internal naming in user-defined steady state files.'])
        problem_dummy=1;
    end
    if strmatch('M_',M.endo_names,'exact') 
        disp(['MODEL_DIAGNOSTICS: using the name M_ for an endogenous variable will typically conflict with the internal naming in user-defined steady state files.'])
        problem_dummy=1;
    end
    if strmatch('M_',M.param_names,'exact')
        disp(['MODEL_DIAGNOSTICS: using the name M_ for a parameter will typically conflict with the internal naming in user-defined steady state files.'])
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

if M.exo_nbr == 0
    oo.exo_steady_state = [] ;
end


info=test_for_deep_parameters_calibration(M);
if info
    problem_dummy=1;
end

% check if ys is steady state
options.debug=true; %locally set debug option to true
if options.logged_steady_state %if steady state was previously logged, undo this
    oo.dr.ys=exp(oo.dr.ys);
    oo.steady_state=exp(oo.steady_state);
    options.logged_steady_state=0;
end
[dr.ys,M.params,check1]=evaluate_steady_state(oo.steady_state,M,options,oo,options.steadystate.nocheck);

% testing for problem
if check1(1)
    problem_dummy=1;
    disp('MODEL_DIAGNOSTICS: The steady state cannot be computed')
    if any(isnan(dr.ys))
        disp(['MODEL_DIAGNOSTICS: Steady state contains NaNs'])
    end
    if any(isinf(dr.ys))
        disp(['MODEL_DIAGNOSTICS: Steady state contains Inf'])
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
if ~isfield(M,'block_structure_stat')
    nb = 1;
else
    nb = length(M.block_structure_stat.block);
end

exo = [oo.exo_steady_state; oo.exo_det_steady_state];
for b=1:nb
    if options.bytecode
        if nb == 1
            [res, jacob] = bytecode(dr.ys, exo, M.params, dr.ys, 1, exo, ...
                                    'evaluate', 'static');
        else
            [res, jacob] = bytecode(dr.ys, exo, M.params, dr.ys, 1, exo, ...
                                    'evaluate', 'static',['block=' ...
                                int2str(b)]);
        end
    else
        [res,jacob]=feval([M.fname '.static'],dr.ys,exo,M.params);
    end
    if any(any(isinf(jacob) | isnan(jacob)))
        problem_dummy=1;
        [infrow,infcol]=find(isinf(jacob) | isnan(jacob));
        fprintf('\nMODEL_DIAGNOSTICS: The Jacobian of the static model contains Inf or NaN. The problem arises from: \n\n')
        display_problematic_vars_Jacobian(infrow,infcol,M,dr.ys,'static','MODEL_DIAGNOSTICS: ')
    end
    if any(any(~isreal(jacob)))
        problem_dummy=1;
        [imagrow,imagcol]=find(abs(imag(jacob))>1e-15);
        fprintf('\nMODEL_DIAGNOSTICS: The Jacobian of the static model contains imaginary parts. The problem arises from: \n\n')
        display_problematic_vars_Jacobian(imagrow,imagcol,M,dr.ys,'static','MODEL_DIAGNOSTICS: ')
    end
    try
        if isoctave || matlab_ver_less_than('9.12') || isempty(options_.jacobian_tolerance)
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
        disp(['MODEL_DIAGNOSTICS:  there is ' num2str(endo_nbr-rank_jacob) ...
              ' colinear relationships between the variables and the equations'])
        if isoctave || matlab_ver_less_than('9.12') || isempty(options_.jacobian_tolerance)
            ncol = null(jacob);
        else
            ncol = null(jacob,options_.jacobian_tolerance); %can sometimes fail
        end
        n_rel = size(ncol,2);
        for i = 1:n_rel
            if n_rel  > 1
                disp(['Relation ' int2str(i)])
            end
            disp('Colinear variables:')
            for j=1:10
                k = find(abs(ncol(:,i)) > 10^-j);
                if max(abs(jacob(:,k)*ncol(k,i))) < 1e-6
                    break
                end
            end
            fprintf('%s\n',endo_names{k})
        end
        if isoctave || matlab_ver_less_than('9.12') || isempty(options_.jacobian_tolerance)
            neq = null(jacob'); %can sometimes fail
        else
            neq = null(jacob',options_.jacobian_tolerance); %can sometimes fail
        end
        n_rel = size(neq,2);
        for i = 1:n_rel
            if n_rel  > 1
                disp(['Relation ' int2str(i)])
            end
            disp('Colinear equations')
            for j=1:10
                k = find(abs(neq(:,i)) > 10^-j);
                if max(abs(jacob(k,:)'*neq(k,i))) < 1e-6
                    break
                end
            end
            disp(k')
        end
    end
end

if singularity_problem
    try
        options_check=options;
        options_check.noprint=1;
        [eigenvalues_] = check(M, options_check, oo);
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
klen = M.maximum_lag + M.maximum_lead + 1;
exo_simul = [repmat(oo.exo_steady_state',klen,1) repmat(oo.exo_det_steady_state',klen,1)];
iyv = M.lead_lag_incidence';
iyv = iyv(:);
iyr0 = find(iyv) ;
it_ = M.maximum_lag + 1;
z = repmat(dr.ys,1,klen);

if ~options.block
    if options.order == 1
        if (options.bytecode)
            [~, loc_dr] = bytecode('dynamic','evaluate', z,exo_simul, ...
                                   M.params, dr.ys, 1);
            jacobia_ = [loc_dr.g1 loc_dr.g1_x loc_dr.g1_xd];
        else
            [~,jacobia_] = feval([M.fname '.dynamic'],z(iyr0),exo_simul, ...
                                 M.params, dr.ys, it_);
        end
    elseif options.order >= 2
        if (options.bytecode)
            [~, loc_dr] = bytecode('dynamic','evaluate', z,exo_simul, ...
                                   M.params, dr.ys, 1);
            jacobia_ = [loc_dr.g1 loc_dr.g1_x];
        else
            [~,jacobia_,hessian1] = feval([M.fname '.dynamic'],z(iyr0),...
                                          exo_simul, ...
                                          M.params, dr.ys, it_);
        end
    end

    if any(any(isinf(jacobia_) | isnan(jacobia_)))
        problem_dummy=1;
        [infrow,infcol]=find(isinf(jacobia_) | isnan(jacobia_));
        fprintf('\nMODEL_DIAGNOSTICS: The Jacobian of the dynamic model contains Inf or NaN. The problem arises from: \n\n')
        display_problematic_vars_Jacobian(infrow,infcol,M,dr.ys,'dynamic','MODEL_DIAGNOSTICS: ')
    end
    if any(any(~isreal(jacobia_)))
        [imagrow,imagcol]=find(abs(imag(jacobia_))>1e-15);
        if ~isempty(imagrow)
            problem_dummy=1;
            fprintf('\nMODEL_DIAGNOSTICS: The Jacobian of the dynamic model contains imaginary parts. The problem arises from: \n\n')
            display_problematic_vars_Jacobian(imagrow,imagcol,M,dr.ys,'dynamic','MODEL_DIAGNOSTICS: ')
        end
    end
    if exist('hessian1','var')
        if any(any(isinf(hessian1) | isnan(hessian1)))
            problem_dummy=1;
            fprintf('\nMODEL_DIAGNOSTICS: The Hessian of the dynamic model contains Inf or NaN.\n')
        end
    end
else
    fprintf('\nMODEL_DIAGNOSTICS: This command currently does not support the block option for checking.\n')
    fprintf('\nMODEL_DIAGNOSTICS: the dynamic model. You may want to disable it for doing model_diagnostics. Skipping this part.\n')
end

if problem_dummy==0
    fprintf('MODEL_DIAGNOSTICS:  No obvious problems with this mod-file were detected.\n')
end
