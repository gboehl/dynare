function [info, oo_, options_, M_] = run(M_, options_, oo_, var_list, params, i_var,W)
% [info, oo_, options_, M_] = osr(M_, options_, oo_, var_list, params, i_var,W)
% Function computing the solution to the optimal simple rule-problem
%
% INPUTS
%   M_          [structure]                 Dynare's model structure
%   oo_         [structure]                 Dynare's results structure
%   options_    [structure]                 Dynare's options structure
%   var_list    [character array]           list of endogenous variables specified
%   params      [character array]           list of parameter to be chosen in
%                                           optimal simple rule
%   i_var       [n_osr_vars by 1 double]    indices of osr-variable in
%                                           specified in optim_weights in declaration order
%   W           [M_.endo_nbr by M_.endo_nbr sparse matrix] Weighting matrix for variance of endogenous variables
%
% OUTPUTS
%   info        [integer]           scalar or vector, error code.
%   oo_         [structure]         Dynare's results structure, containing subfield   
%       osr_res:    results structure containing:
%                   - objective_function [scalar double]   value of the objective
%                                                           function at the optimum
%                   - optim_params       [structure]       parameter values at the optimum
%   options_    [structure]         Dynare's options structure
%   M_          [structure]         Dynare's model structure
%
% SPECIAL REQUIREMENTS
%   none.

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

options_.order = 1;

if isempty(options_.qz_criterium)
    options_.qz_criterium = 1+1e-6;
end

oo_=make_ex_(M_,options_,oo_);

np = size(params,1);
i_params = zeros(np,1);
for i=1:np
    str = deblank(params(i,:));
    i_params(i) = strmatch(str{:}, M_.param_names, 'exact');
end

if ~options_.noprint
    fprintf('\nOPTIMAL SIMPLE RULE\n')
end

osr_res.error_indicator = 1; %initialize indicator

if M_.exo_nbr == 0
    oo_.exo_steady_state = [] ;
end

if ~M_.lead_lag_incidence(M_.maximum_lag+1,:) > 0
    error ('OSR: Error in model specification: some variables don''t appear as current') ;
end

if M_.maximum_lead == 0
    error ('OSR: Backward or static model: no point in using OSR') ;
end

if any(any(isinf(W)))
    error ('OSR: At least one of the optim_weights is infinite.') ;
end

if any(isnan(M_.params(i_params)))
    error ('OSR: At least one of the initial parameter values for osr_params is NaN') ;
end

%restore backward compatibility with maxit and tolf
if isfield(options_.osr,'maxit') || isfield(options_.osr,'tolf')
    warning('OSR: The use of maxit and tolf is now deprecated. Please use the optim-option instead.')
    if options_.osr.opt_algo~=4
        error ('OSR: The maxit and tolf options are not supported when not using the default opt_algo=4. Use the optim-option instead.') ;
    else
        if isfield(options_.osr,'maxit')
            if ~isempty(options_.optim_opt)
                options_.optim_opt=[options_.optim_opt,','];
            end
            options_.optim_opt=[options_.optim_opt,'''MaxIter'',',num2str(options_.osr.maxit),''];
        end
        if isfield(options_.osr,'tolf')
            if ~isempty(options_.optim_opt)
                options_.optim_opt=[options_.optim_opt,','];
            end
            options_.optim_opt=[options_.optim_opt,'''TolFun'',',num2str(options_.osr.tolf),''];
        end
    end
end

oo_.dr = set_state_space(oo_.dr,M_);
par_0 = M_.params(i_params);
inv_order_var = oo_.dr.inv_order_var;

%extract unique entries of covariance
i_var=unique(i_var);
%% do initial checks
[loss,info]=osr.objective(par_0,M_,oo_,options_,i_params,inv_order_var(i_var),W(i_var,i_var));
if info~=0
    print_info(info, options_.noprint, options_);
else
    if ~options_.noprint
        fprintf('\nOSR: Initial value of the objective function: %g \n\n',loss);
    end
end
if ~options_.noprint && isinf(loss)
    fprintf('\nOSR: The initial value of the objective function is infinite.\n');
    fprintf('\nOSR: Check whether the unconditional variance of a target variable is infinite\n');
    fprintf('\nOSR: due to the presence of a unit root.\n');
    error('OSR: Initial likelihood is infinite')
end

if isequal(options_.osr.opt_algo,5)
    error('OSR: OSR does not support opt_algo=5.')
elseif isequal(options_.osr.opt_algo,6)
    error('OSR: OSR does not support opt_algo=6.')
elseif isequal(options_.osr.opt_algo,10)
    error('OSR: OSR does not support opt_algo=10.')
elseif isequal(options_.osr.opt_algo,11)
    error('OSR: OSR does not support opt_algo=11.')
else
    if ~isempty(M_.osr.param_bounds) && ~(ismember(options_.osr.opt_algo,[1,2,5,9]) || ischar(options_.osr.opt_algo))
        error('OSR: OSR with bounds on parameters requires a constrained optimizer, i.e. opt_algo= 1,2 or 9.')
    end
    %%do actual optimization
    [p, f] = dynare_minimize_objective(str2func('osr.objective'),par_0,options_.osr.opt_algo,options_,M_.osr.param_bounds,M_.param_names(i_params),[],[], M_,oo_,options_,i_params,...
                                       inv_order_var(i_var),W(i_var,i_var));
end

osr_res.objective_function = f;
M_.params(i_params)=p; %make sure optimal parameters are set (and not the last draw used in csminwel)
for i=1:length(i_params)
    osr_res.optim_params.(M_.param_names{i_params(i)}) = p(i);
end

if ~options_.noprint
    my_title='OPTIMAL VALUE OF THE PARAMETERS';
    labels = M_.param_names(i_params);
    headers = {'Variables'; 'Value'};
    lh = cellofchararraymaxlength(labels)+2;
    dyntable(options_, my_title, headers, labels, p, lh, 10, 6);
    if options_.TeX
        labels = M_.param_names_tex(i_params);       
        lh = cellofchararraymaxlength(labels)+2;
        dyn_latex_table(M_, options_, my_title, 'osr', headers, labels, p, lh, 10, 6);
    end
    fprintf('\nObjective function : %16.6g\n\n',f);
end

[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list);

if ~info(1)
    osr_res.error_indicator=0;
end
oo_.osr = osr_res;