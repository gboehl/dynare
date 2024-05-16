function [endogenousvariables, success, err, iter] = sim1(endogenousvariables, exogenousvariables, steadystate, M_, options_)
% [endogenousvariables, success, err, iter] = sim1(endogenousvariables, exogenousvariables, steadystate, M_, options_)
% Performs deterministic simulations with lead or lag of one period, using
% a basic Newton solver on sparse matrices.
% Uses perfect_foresight_problem DLL to construct the stacked problem.
%
% INPUTS
%   - endogenousvariables [double] N*(T+M_.maximum_lag+M_.maximum_lead) array, paths for the endogenous variables (initial condition + initial guess + terminal condition).
%   - exogenousvariables  [double] T*M array, paths for the exogenous variables.
%   - steadystate         [double] N*1 array, steady state for the endogenous variables.
%   - M_                  [struct] contains a description of the model.
%   - options_            [struct] contains various options.
% OUTPUTS
%   - endogenousvariables [double] N*(T+M_.maximum_lag+M_.maximum_lead) array, paths for the endogenous variables (solution of the perfect foresight model).
%   - success             [logical] Whether a solution was found
%   - err                 [double] ∞-norm of the residual
%   - iter                [integer] Number of iterations

% Copyright © 1996-2023 Dynare Team
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

verbose = options_.verbosity && ~options_.noprint;

ny = M_.endo_nbr;
periods = options_.periods;
vperiods = periods*ones(1,options_.simul.maxit);

if M_.maximum_lag > 0
    y0 = endogenousvariables(:, M_.maximum_lag);
else
    y0 = NaN(ny, 1);
end

if M_.maximum_lead > 0
    yT = endogenousvariables(:, M_.maximum_lag+periods+1);
else
    yT = NaN(ny, 1);
end

y = reshape(endogenousvariables(:, M_.maximum_lag+(1:periods)), ny*periods, 1);

stop = false;

if verbose
    skipline()
    printline(56)
    disp('MODEL SIMULATION:')
    skipline()
end

h1 = clock;

for iter = 1:options_.simul.maxit
    h2 = clock;

    [res, A] = perfect_foresight_problem(y, y0, yT, exogenousvariables, M_.params, steadystate, periods, M_, options_);
    % A is the stacked Jacobian with period x equations alongs the rows and
    % periods times variables (in declaration order) along the columns
    if options_.debug && iter==1
        [row,col]=find(A);
        row=setdiff(1:periods*ny,row);
        column=setdiff(1:periods*ny,col);
        if ~isempty(row) || ~isempty(column)
            fprintf('The stacked Jacobian is singular. The problem derives from:\n')
            if ~isempty(row)
                time_period=ceil(row/ny);
                equation=row-ny*(time_period-1);
                for eq_iter=1:length(equation)
                    fprintf('The derivative of equation %d at time %d is zero for all variables\n',equation(eq_iter),time_period(eq_iter));
                end
            end
            if ~isempty(column)
                time_period=ceil(column/ny);
                variable=column-ny*(time_period-1);
                for eq_iter=1:length(variable)
                    fprintf('The derivative with respect to variable %d at time %d is zero for all equations\n',variable(eq_iter),time_period(eq_iter));
                end
            end            
        end
        check_Jacobian_for_singularity(full(A),M_.endo_names,options_);
    end
    if options_.endogenous_terminal_period && iter > 1
        for it = 1:periods
            if max(abs(res((it-1)*ny+(1:ny)))) < options_.dynatol.f/1e7
                if it < periods
                    res = res(1:(it*ny));
                    A = A(1:(it*ny), 1:(it*ny));
                    yT = y(it*ny+(1:ny));
                    endogenousvariables(:, M_.maximum_lag+((it+1):periods)) = reshape(y(it*ny+(1:(ny*(periods-it)))), ny, periods-it);
                    y = y(1:(it*ny));
                    periods = it;
                end
                break
            end
        end
        vperiods(iter) = periods;
    end

    err = max(abs(res));
    if options_.debug
        fprintf('\nLargest absolute residual at iteration %d: %10.3f\n',iter,err);
        if any(isnan(res)) || any(isinf(res)) || any(any(isnan(endogenousvariables))) || any(any(isinf(endogenousvariables)))
            fprintf('\nWARNING: NaN or Inf detected in the residuals or endogenous variables.\n');
        end
        skipline()
    end
    if verbose
        fprintf('Iter: %d,\t err. = %g,\t time = %g\n', iter, err, etime(clock,h2));
    end
    if err < options_.dynatol.f
        stop = true;
        break
    end
    if options_.simul.robust_lin_solve
        dy = -lin_solve_robust(A, res, verbose, options_);
    else
        dy = -lin_solve(A, res, verbose);
    end
    if any(isnan(dy)) || any(isinf(dy))
        if verbose
            display_critical_variables(reshape(dy,[ny periods])', M_, options_.noprint);
        end
    end
    y = y + dy;
end

endogenousvariables(:, M_.maximum_lag+(1:periods)) = reshape(y, ny, periods);

if options_.endogenous_terminal_period
    periods = options_.periods;
    err = evaluate_max_dynamic_residual(str2func([M_.fname,'.dynamic']), endogenousvariables, exogenousvariables, M_.params, steadystate, periods, ny, M_.maximum_endo_lag, M_.lead_lag_incidence);
end

if stop
    % initial or terminal observations may contain
    % harmless NaN or Inf. We test only values computed above
    if any(any(isnan(y))) || any(any(isinf(y)))
        success = false; % NaN or Inf occurred
        if verbose
            skipline()
            fprintf('Total time of simulation: %g.\n', etime(clock,h1))
            disp('Simulation terminated with NaN or Inf in the residuals or endogenous variables.')
            display_critical_variables(reshape(dy,[ny periods])', M_, options_.noprint);
            disp('There is most likely something wrong with your model. Try model_diagnostics or another simulation method.')
            printline(105)
        end
    else
        if verbose
            skipline();
            fprintf('Total time of simulation: %g.\n', etime(clock,h1))
            printline(56)
        end
        success = true; % Convergency obtained.
    end
elseif ~stop
    if verbose
        skipline();
        fprintf('Total time of simulation: %g.\n', etime(clock,h1))
        disp('Maximum number of iterations is reached (modify option maxit).')
        printline(62)
    end
    success = false; % more iterations are needed.
end

if verbose
    skipline();
end

function x = lin_solve(A, b, verbose)
if norm(b) < sqrt(eps) % then x = 0 is a solution
    x = 0;
    return
end

x = A\b;
x(~isfinite(x)) = 0;
relres = norm(b - A*x) / norm(b);
if relres > 1e-6 && verbose
    fprintf('WARNING : Failed to find a solution to the linear system.\n');
end

function [ x, flag, relres ] = lin_solve_robust(A, b ,verbose, options_)
if norm(b) < sqrt(eps) % then x = 0 is a solution
    x = 0;
    flag = 0;
    relres = 0;
    return
end

x = A\b;
x(~isfinite(x)) = 0;
[ x, flag, relres ] = bicgstab(A, b, [], [], [], [], x); % returns immediately if x is a solution
if flag == 0
    return
end

if ~options_.noprint
    disp( relres );
end

if verbose
    fprintf('Initial bicgstab failed, trying alternative start point.\n');
end
old_x = x;
old_relres = relres;
[ x, flag, relres ] = bicgstab(A, b);
if flag == 0
    return
end

if verbose
    fprintf('Alternative start point also failed with bicgstab, trying gmres.\n');
end
if old_relres < relres
    x = old_x;
end
[ x, flag, relres ] = gmres(A, b, [], [], [], [], [], x);
if flag == 0
    return
end

if verbose
    fprintf('Initial gmres failed, trying alternative start point.\n');
end
old_x = x;
old_relres = relres;
[ x, flag, relres ] = gmres(A, b);
if flag == 0
    return
end

if verbose
    fprintf('Alternative start point also failed with gmres, using the (SLOW) Moore-Penrose Pseudo-Inverse.\n');
end
if old_relres < relres
    x = old_x;
    relres = old_relres;
end
old_x = x;
old_relres = relres;
x = pinv(full(A)) * b;
relres = norm(b - A*x) / norm(b);
if old_relres < relres
    x = old_x;
    relres = old_relres;
end
flag = relres > 1e-6;
if flag ~= 0 && verbose
    fprintf('WARNING : Failed to find a solution to the linear system\n');
end

function display_critical_variables(dyy, M_, noprint)

if noprint
    return
end

if any(isnan(dyy))
    indx = find(any(isnan(dyy)));
    endo_names= M_.endo_names(indx);
    disp('Last iteration provided NaN for the following variables:')
    fprintf('%s, ', endo_names{:}),
    fprintf('\n'),
end
if any(isinf(dyy))
    indx = find(any(isinf(dyy)));
    endo_names = M_.endo_names(indx);
    disp('Last iteration diverged (Inf) for the following variables:')
    fprintf('%s, ', endo_names{:}),
    fprintf('\n'),
end


function check_Jacobian_for_singularity(jacob,endo_names,options_)
n_vars_jacob=size(jacob,2);
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
    disp(['sim1:  The Jacobian of the dynamic model is ' ...
        'singular'])
    disp(['sim1:  there is ' num2str(n_vars_jacob-rank_jacob) ...
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
        fprintf('%s\n',endo_names{mod(k-1,length(endo_names))+1})
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
        equation=mod(k-1,length(endo_names))+1;
        period=ceil(k/length(endo_names));
        for ii=1:length(equation)
            fprintf('Equation %5u, period %5u\n',equation(ii),period(ii))
        end
    end
end