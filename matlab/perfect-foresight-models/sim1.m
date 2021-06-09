function [endogenousvariables, info] = sim1(endogenousvariables, exogenousvariables, steadystate, M, options)
% Performs deterministic simulations with lead or lag of one period, using
% a basic Newton solver on sparse matrices.
% Uses perfect_foresight_problem DLL to construct the stacked problem.
%
% INPUTS
%   - endogenousvariables [double] N*(T+M.maximum_lag+M.maximum_lead) array, paths for the endogenous variables (initial condition + initial guess + terminal condition).
%   - exogenousvariables  [double] T*M array, paths for the exogenous variables.
%   - steadystate         [double] N*1 array, steady state for the endogenous variables.
%   - M                   [struct] contains a description of the model.
%   - options             [struct] contains various options.
% OUTPUTS
%   - endogenousvariables [double] N*(T+M.maximum_lag+M.maximum_lead) array, paths for the endogenous variables (solution of the perfect foresight model).
%   - info                [struct] contains informations about the results.

% Copyright (C) 1996-2021 Dynare Team
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

verbose = options.verbosity && ~options.noprint;

ny = M.endo_nbr;
periods = options.periods;
vperiods = periods*ones(1,options.simul.maxit);

if M.maximum_lag > 0
    y0 = endogenousvariables(:, M.maximum_lag);
else
    y0 = NaN(ny, 1);
end

if M.maximum_lead > 0
    yT = endogenousvariables(:, M.maximum_lag+periods+1);
else
    yT = NaN(ny, 1);
end

y = reshape(endogenousvariables(:, M.maximum_lag+(1:periods)), ny*periods, 1);

stop = false;

if verbose
    skipline()
    printline(56)
    disp('MODEL SIMULATION:')
    skipline()
end

h1 = clock;

for iter = 1:options.simul.maxit
    h2 = clock;

    [res, A] = perfect_foresight_problem(y, y0, yT, exogenousvariables, M.params, steadystate, periods, M, options);
    % A is the stacked Jacobian with period x equations alongs the rows and
    % periods times variables (in declaration order) along the columns
    if options.debug && iter==1
        row=find(all(A==0,2));
        column=find(all(A==0,1));
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
    end
    if options.endogenous_terminal_period && iter > 1
        for it = 1:periods
            if max(abs(res((it-1)*ny+(1:ny)))) < options.dynatol.f/1e7
                if it < periods
                    res = res(1:(it*ny));
                    A = A(1:(it*ny), 1:(it*ny));
                    yT = y(it*ny+(1:ny));
                    endogenousvariables(:, M.maximum_lag+((it+1):periods)) = reshape(y(it*ny+(1:(ny*(periods-it)))), ny, periods-it);
                    y = y(1:(it*ny));
                    periods = it;
                end
                break
            end
        end
        vperiods(iter) = periods;
    end

    err = max(abs(res));
    if options.debug
        fprintf('\nLargest absolute residual at iteration %d: %10.3f\n',iter,err);
        if any(isnan(res)) || any(isinf(res)) || any(any(isnan(endogenousvariables))) || any(any(isinf(endogenousvariables)))
            fprintf('\nWARNING: NaN or Inf detected in the residuals or endogenous variables.\n');
        end
        skipline()
    end
    if verbose
        fprintf('Iter: %d,\t err. = %g,\t time = %g\n', iter, err, etime(clock,h2));
    end
    if err < options.dynatol.f
        stop = true;
        break
    end
    if options.simul.robust_lin_solve
        dy = -lin_solve_robust(A, res, verbose, options);
    else
        dy = -lin_solve(A, res, verbose);
    end
    if any(isnan(dy)) || any(isinf(dy))
        if verbose
            display_critical_variables(reshape(dy,[ny periods])', M, options.noprint);
        end
    end
    y = y + dy;
end

endogenousvariables(:, M.maximum_lag+(1:periods)) = reshape(y, ny, periods);

if options.endogenous_terminal_period
    periods = options.periods;
    err = evaluate_max_dynamic_residual(str2func([M.fname,'.dynamic']), endogenousvariables, exogenousvariables, M.params, steadystate, periods, ny, M.maximum_endo_lag, M.lead_lag_incidence);
end

if stop
    % initial or terminal observations may contain
    % harmless NaN or Inf. We test only values computed above
    if any(any(isnan(y))) || any(any(isinf(y)))
        info.status = false;% NaN or Inf occurred
        info.error = err;
        info.iterations = iter;
        info.periods = vperiods(1:iter);
        if verbose
            skipline()
            fprintf('Total time of simulation: %g.\n', etime(clock,h1))
            disp('Simulation terminated with NaN or Inf in the residuals or endogenous variables.')
            display_critical_variables(reshape(dy,[ny periods])', M, options.noprint);
            disp('There is most likely something wrong with your model. Try model_diagnostics or another simulation method.')
            printline(105)
        end
    else
        if verbose
            skipline();
            fprintf('Total time of simulation: %g.\n', etime(clock,h1))
            printline(56)
        end
        info.status = true;% Convergency obtained.
        info.error = err;
        info.iterations = iter;
        info.periods = vperiods(1:iter);
    end
elseif ~stop
    if verbose
        skipline();
        fprintf('Total time of simulation: %g.\n', etime(clock,h1))
        disp('Maximum number of iterations is reached (modify option maxit).')
        printline(62)
    end
    info.status = false;% more iterations are needed.
    info.error = err;
    info.periods = vperiods(1:iter);
    info.iterations = options.simul.maxit;
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

function [ x, flag, relres ] = lin_solve_robust(A, b ,verbose, options)
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

if ~options.noprint
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

function display_critical_variables(dyy, M, noprint)

if noprint
    return
end

if any(isnan(dyy))
    indx = find(any(isnan(dyy)));
    endo_names= M.endo_names(indx);
    disp('Last iteration provided NaN for the following variables:')
    fprintf('%s, ', endo_names{:}),
    fprintf('\n'),
end
if any(isinf(dyy))
    indx = find(any(isinf(dyy)));
    endo_names = M.endo_names(indx);
    disp('Last iteration diverged (Inf) for the following variables:')
    fprintf('%s, ', endo_names{:}),
    fprintf('\n'),
end
