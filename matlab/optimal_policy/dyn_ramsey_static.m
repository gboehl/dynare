function [steady_state, params, check] = dyn_ramsey_static(ys_init, exo_ss, M_, options_)
% [steady_state, params, check] = dyn_ramsey_static(ys_init, exo_ss, M_, options_)
% Computes the steady state for optimal policy
%
% When there is no steady state file, relies on the fact that Lagrange
% multipliers appear linearly in the system to be solved. Instead of directly
% solving for the Lagrange multipliers along with the other variables, the
% algorithms reduces the size of the problem by always computing the value of
% the multipliers that minimizes the residuals, given the other variables
% (using a minimum norm solution, easy to compute because of the linearity
% property).
%
% INPUTS
%    ys_init:       vector of endogenous variables or instruments
%    exo_ss         vector of exogenous steady state (incl. deterministic exogenous)
%    M_:             Dynare model structure
%    options:       Dynare options structure
%
% OUTPUTS
%    steady_state:  steady state value
%    params:        parameters at steady state, potentially updated by
%                   steady_state file
%    check:         error indicator, 0 if everything is OK
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2003-2024 Dynare Team
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


params = M_.params;
check = 0;
options_.steadystate.nocheck = 1; %locally disable checking because Lagrange multipliers are not accounted for in evaluate_steady_state_file
nl_func = @(x) dyn_ramsey_static_1(x,exo_ss,ys_init,M_,options_);

if ~options_.steadystate_flag && check_static_model(ys_init,exo_ss,M_,options_)
    steady_state = ys_init;
    return
elseif options_.steadystate_flag
    k_inst = [];
    inst_nbr = size(options_.instruments,1);
    for i = 1:inst_nbr
        k_inst = [k_inst; strmatch(options_.instruments{i}, M_.endo_names, 'exact')];
    end
    if inst_nbr == 1
        %solve for instrument, using univariate solver, starting at initial value for instrument
        [inst_val, info1]= csolve(nl_func,ys_init(k_inst),'',options_.solve_tolf,options_.ramsey.maxit);
        if info1==1 || info1==3
            check=81;
        end
        if info1==4
            check=87;
        end
    else
        %solve for instrument, using multivariate solver, starting at
        %initial value for instruments
        o_jacobian_flag = options_.jacobian_flag;
        options_.jacobian_flag = false;
        [inst_val, errorflag] = dynare_solve(nl_func, ys_init(k_inst), options_.ramsey.maxit, options_.solve_tolf, options_.solve_tolx, options_);
        options_.jacobian_flag = o_jacobian_flag;
        if errorflag
            check=81;
        end
    end
    ys_init(k_inst) = inst_val;
    [~,params] = evaluate_steady_state_file(ys_init,exo_ss,M_,options_,~options_.steadystate.nocheck); %run steady state file again to update parameters
    [~,~,steady_state] = nl_func(inst_val); %compute and return steady state
else
    xx = ys_init(1:M_.orig_endo_nbr);
    o_jacobian_flag = options_.jacobian_flag;
    options_.jacobian_flag = false;
    [xx, errorflag] = dynare_solve(nl_func, xx, options_.ramsey.maxit, options_.solve_tolf, options_.solve_tolx, options_);
    options_.jacobian_flag = o_jacobian_flag;
    if errorflag
        check=81;
    end
    [~,~,steady_state] = nl_func(xx);
end


function [resids,rJ,steady_state] = dyn_ramsey_static_1(x,exo_ss,ys_init,M_,options_)
resids = [];
rJ = [];
mult = [];

inst_nbr = M_.ramsey_orig_endo_nbr - M_.ramsey_orig_eq_nbr;

if options_.steadystate_flag
    k_inst = [];
    for i = 1:size(options_.instruments,1)
        k_inst = [k_inst; strmatch(options_.instruments{i}, M_.endo_names, 'exact')];
    end
    ys_init(k_inst) = x; %set instrument, the only value required for steady state computation, to current value
    [x,M_.params,check] = evaluate_steady_state_file(ys_init,... %returned x now has size endo_nbr as opposed to input size of n_instruments
                                                    exo_ss, ...
                                                    M_,options_,~options_.steadystate.nocheck);
    if any(imag(x(1:M_.orig_endo_nbr))) %return with penalty
        resids=ones(inst_nbr,1)+sum(abs(imag(x(1:M_.orig_endo_nbr)))); %return with penalty
        steady_state=NaN(M_.endo_nbr,1);
        return
    end
    if check(1) %return 
        resids=ones(inst_nbr,1)+sum(abs(x(1:M_.orig_endo_nbr))); %return with penalty
        steady_state=NaN(M_.endo_nbr,1);
        return        
    end
end

xx = zeros(M_.endo_nbr,1); %initialize steady state vector
xx(1:M_.orig_endo_nbr) = x(1:M_.orig_endo_nbr); %set values of original endogenous variables based on steady state file or initial value

% Determine whether other auxiliary variables will need to be updated
if any([M_.aux_vars.type] ~= 6) %auxiliary variables other than multipliers
    needs_set_auxiliary_variables = true;
    fh = str2func([M_.fname '.set_auxiliary_variables']);
    s_a_v_func = @(z) fh(z, exo_ss, M_.params);
    xx = s_a_v_func(xx);
else
    needs_set_auxiliary_variables = false;
end

% Compute the value of the Lagrange multipliers that minimizes the norm of the
% residuals, given the other endogenous
if options_.bytecode
    res = bytecode('static', M_, options_, xx, exo_ss, M_.params, 'evaluate');
else
    res = feval([M_.fname '.sparse.static_resid'], xx, exo_ss, M_.params);
end
A = feval([M_.fname '.ramsey_multipliers_static_g1'], xx, exo_ss, M_.params, M_.ramsey_multipliers_static_g1_sparse_rowval, M_.ramsey_multipliers_static_g1_sparse_colval, M_.ramsey_multipliers_static_g1_sparse_colptr);
y = res(1:M_.ramsey_orig_endo_nbr);
mult = -A\y;

resids1 = y+A*mult;
if inst_nbr == 1
    r1 = sqrt(resids1'*resids1);
else
    [~,r] = qr([A y]');
    k = size(A,1)+(1-inst_nbr:0);
    r1 = r(end,k)';
end
if options_.steadystate_flag
    resids = r1;
else
    resids = [res(M_.ramsey_orig_endo_nbr+(1:M_.orig_endo_nbr-inst_nbr)); r1];
end
if needs_set_auxiliary_variables
    steady_state = s_a_v_func([xx(1:M_.ramsey_orig_endo_nbr); mult]);
else
    steady_state = [xx(1:M_.ramsey_orig_endo_nbr); mult];
end

function result = check_static_model(ys,exo_ss,M_,options_)
result = false;
if (options_.bytecode)
    res = bytecode('static', M_, options_, ys, exo_ss, M_.params, 'evaluate');
else
    res = feval([M_.fname '.sparse.static_resid'], ys, exo_ss, M_.params);
end
if norm(res) < options_.solve_tolf
    result = true;
end
