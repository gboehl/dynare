function [x, errorflag, fvec, fjac, errorcode] = dynare_solve(f, x, maxit, tolf, tolx, options, varargin)

% Solves a nonlinear system of equations, f(x) = 0 with n unknowns
% and n equations.
%
% INPUTS
% - f            [char, fhandle]    function to be solved
% - x            [double]           n×1 vector, initial guess.
% - options      [struct]           Dynare options, aka options_.
% - varargin                        list of additional arguments to be passed to func.
%
% OUTPUTS
% - x            [double]           n×1 vector, solution.
% - errorflag    [logical]          scalar, true iff the model can not be solved.
% - fvec         [double]           n×1 vector, function value at x (f(x), used for debugging when errorflag is true).
% - fjac         [double]           n×n matrix, Jacobian value at x (J(x), used for debugging when errorflag is true).
% - errorcode    [integer]          scalar.
%
% REMARKS
% Interpretation of the error code depends on the algorithm, except if value of errorcode is
%
%        -10  -> System of equation ill-behaved at the initial guess (Inf, Nans or complex numbers).
%        -11  -> Initial guess is a solution of the system of equations.

% Copyright © 2001-2022 Dynare Team
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

jacobian_flag = options.jacobian_flag; % true iff Jacobian is returned by f routine (as a second output argument).

errorflag = false; % Let's be optimistic!
nn = size(x,1);

% Keep a copy of the initial guess.
x0 = x;

% Get status of the initial guess (default values?)
if any(x)
    % The current initial guess is not the default for all the variables.
    idx = find(x);      % Indices of the variables with default initial guess values.
    in0 = length(idx);
else
    % The current initial guess is the default for all the variables.
    idx = transpose(1:nn);
    in0 = nn;
end

% Get first element of varargin if solve_algo ∈ {12,14} and rename varargin.
if ismember(options.solve_algo, [12, 14])
    isloggedlhs = varargin{1};
    isauxdiffloggedrhs = varargin{2};
    endo_names = varargin{3};
    lhs = varargin{4};
    args = varargin(5:end);
else
    args = varargin;
end

% checking initial values
% TODO We should have an option to deactivate the randomization.
if jacobian_flag
    [fvec, fjac] = feval(f, x, args{:});
    wrong_initial_guess_flag = false;
    if ~all(isfinite(fvec)) || any(isinf(fjac(:))) || any(isnan((fjac(:)))) || any(~isreal(fvec)) || any(~isreal(fjac(:)))
        if ~ismember(options.solve_algo,[10,11]) && max(abs(fvec))< tolf
            % return if initial value solves the problem except if a mixed complementarity problem is to be solved (complementarity conditions may not be satisfied)
            errorcode = -11;
            return;
        end
        disp_verbose('Randomize initial guess...', options.verbosity)
        % Let's try random numbers for the variables initialized with the default value.
        wrong_initial_guess_flag = true;
        % First try with positive numbers.
        tentative_number = 0;
        while wrong_initial_guess_flag && tentative_number<=in0*10
            tentative_number = tentative_number+1;
            x(idx) = rand(in0, 1)*10;
            [fvec, fjac] = feval(f, x, args{:});
            wrong_initial_guess_flag = ~all(isfinite(fvec)) || any(isinf(fjac(:))) || any(isnan((fjac(:))));
        end
        % If all previous attempts failed, try with real numbers.
        tentative_number = 0;
        while wrong_initial_guess_flag && tentative_number<=in0*10
            tentative_number = tentative_number+1;
            x(idx) = randn(in0, 1)*10;
            [fvec, fjac] = feval(f, x, args{:});
            wrong_initial_guess_flag = ~all(isfinite(fvec)) || any(isinf(fjac(:))) || any(isnan((fjac(:))));
        end
        % Last tentative, ff all previous attempts failed, try with negative numbers.
        tentative_number = 0;
        while wrong_initial_guess_flag && tentative_number<=in0*10
            tentative_number = tentative_number+1;
            x(idx) = -rand(in0, 1)*10;
            [fvec, fjac] = feval(f, x, args{:});
            wrong_initial_guess_flag = ~all(isfinite(fvec)) || any(isinf(fjac(:))) || any(isnan((fjac(:))));
        end
    end
else
    fvec = feval(f, x, args{:});
    fjac = zeros(nn, nn);
    if ~ismember(options.solve_algo,[10,11]) && max(abs(fvec)) < tolf
        % return if initial value solves the problem except if a mixed complementarity problem is to be solved (complementarity conditions may not be satisfied)
        errorcode = -11;
        return;
    end
    wrong_initial_guess_flag = false;
    if ~all(isfinite(fvec))
        % Let's try random numbers for the variables initialized with the default value.
        wrong_initial_guess_flag = true;
        % First try with positive numbers.
        tentative_number = 0;
        while wrong_initial_guess_flag && tentative_number<=in0*10
            tentative_number = tentative_number+1;
            x(idx) = rand(in0, 1)*10;
            fvec = feval(f, x, args{:});
            wrong_initial_guess_flag = ~all(isfinite(fvec));
        end
        % If all previous attempts failed, try with real numbers.
        tentative_number = 0;
        while wrong_initial_guess_flag && tentative_number<=in0*10
            tentative_number = tentative_number+1;
            x(idx) = randn(in0, 1)*10;
            fvec = feval(f, x, args{:});
            wrong_initial_guess_flag = ~all(isfinite(fvec));
        end
        % Last tentative, ff all previous attempts failed, try with negative numbers.
        tentative_number = 0;
        while wrong_initial_guess_flag && tentative_number<=in0*10
            tentative_number = tentative_number+1;
            x(idx) = -rand(in0, 1)*10;
            fvec = feval(f, x, args{:});
            wrong_initial_guess_flag = ~all(isfinite(fvec));
        end
    end
end

% Exit with error if no initial guess has been found.
if wrong_initial_guess_flag
    errorcode = -10;
    errorflag = true;
    x = x0;
    return
end

if options.solve_algo == 0
    if ~isoctave
        if ~user_has_matlab_license('optimization_toolbox')
            error('You can''t use solve_algo=0 since you don''t have MATLAB''s Optimization Toolbox')
        end
    end
    options4fsolve = optimset('fsolve');
    options4fsolve.MaxFunEvals = 50000;
    options4fsolve.MaxIter = maxit;
    options4fsolve.TolFun = tolf;
    if options.debug==1
        options4fsolve.Display = 'final';
    else
        options4fsolve.Display = 'off';
    end
    if jacobian_flag
        options4fsolve.Jacobian = 'on';
    else
        options4fsolve.Jacobian = 'off';
    end
    if ~isoctave
        [x, ~, errorcode] = fsolve(f, x, options4fsolve, args{:});
    else
        % Under Octave, use a wrapper, since fsolve() does not have a 4th arg
        if ischar(f)
            f2 = str2func(f);
        else
            f2 = f;
        end
        f = @(x) f2(x, args{:});
        [x, ~, errorcode] = fsolve(f, x, options4fsolve);
    end
    if errorcode==1
        errorflag = false;
    elseif errorcode>1
        if ~isoctave
            if ischar(f)
                f2 = str2func(f);
            else
                f2 = f;
            end
            f = @(x) f2(x, args{:});
        end
        if max(abs(fvec)) > tolf
            errorflag = true;
        else
            errorflag = false;
        end
    else
        errorflag = true;
    end
    [fvec, fjac] = feval(f, x, args{:});
elseif options.solve_algo==1
    [x, errorflag, errorcode] = solve1(f, x, 1:nn, 1:nn, jacobian_flag, options.gstep, tolf, tolx, maxit, [], options.debug, args{:});
    [fvec, fjac] = feval(f, x, args{:});
elseif options.solve_algo==9
    [x, errorflag, errorcode] = trust_region(f, x, 1:nn, 1:nn, jacobian_flag, options.gstep, tolf, tolx, maxit, options.trust_region_initial_step_bound_factor, options.debug, args{:});
    [fvec, fjac] = feval(f, x, args{:});
elseif ismember(options.solve_algo, [2, 12, 4])
    if ismember(options.solve_algo, [2, 12])
        solver = @solve1;
    else
        solver = @trust_region;
    end
    specializedunivariateblocks = options.solve_algo == 12;
    if ~jacobian_flag
        fjac = zeros(nn,nn) ;
        dh = max(abs(x), options.gstep(1)*ones(nn,1))*eps^(1/3);
        for j = 1:nn
            xdh = x ;
            xdh(j) = xdh(j)+dh(j) ;
            fjac(:,j) = (feval(f, xdh, args{:})-fvec)./dh(j) ;
        end
    end
    [j1,j2,r,s] = dmperm(fjac);
    JAC = abs(fjac(j1,j2))>0;
    if options.debug
        disp(['DYNARE_SOLVE (solve_algo=2|4|12): number of blocks = ' num2str(length(r)-1)]);
    end
    l = 0;
    fre = false;
    for i=length(r)-1:-1:1
        blocklength = r(i+1)-r(i);
        if options.debug
            dprintf('DYNARE_SOLVE (solve_algo=2|4|12): solving block %u of size %u.', i, blocklength);
        end
        j = r(i):r(i+1)-1;
        if specializedunivariateblocks
            if options.debug
                dprintf('DYNARE_SOLVE (solve_algo=2|4|12): solving block %u by evaluating RHS.', i);
            end
            if isequal(blocklength, 1)
                if i<length(r)-1
                    if fre || any(JAC(r(i), s(i)+(1:l)))
                        % Reevaluation of the residuals is required because the current RHS depends on
                        % variables that potentially have been updated previously.
                        z = feval(f, x, args{:});
                        l = 0;
                        fre = false;
                    end
                else
                    % First iteration requires the evaluation of the residuals.
                    z = feval(f, x, args{:});
                end
                l = l+1;
                if isequal(lhs{j1(j)}, endo_names{j2(j)}) || isequal(lhs{j1(j)}, sprintf('log(%s)', endo_names{j2(j)}))
                    if isloggedlhs(j1(j))
                        x(j2(j)) = exp(log(x(j2(j)))-z(j1(j)));
                    else
                        x(j2(j)) = x(j2(j))-z(j1(j));
                    end
                else
                    if options.debug
                        dprintf('LHS variable is not determined by RHS expression (%u).', j1(j))
                        dprintf('%s -> %s', lhs{j1(j)}, endo_names{j2(j)})
                    end
                    if ~isempty(regexp(lhs{j1(j)}, '\<AUX_DIFF_(\d*)\>', 'once'))
                        if isauxdiffloggedrhs(j1(j))
                            x(j2(j)) = exp(log(x(j2(j)))+z(j1(j)));
                        else
                            x(j2(j)) = x(j2(j))+z(j1(j));
                        end
                    else
                        error('Algorithm solve_algo=%u cannot be used with this nonlinear problem.', options.solve_algo)
                    end
                end
                continue
            end
        else
            if options.debug
                dprintf('DYNARE_SOLVE (solve_algo=2|4|12): solving block %u with trust_region routine.', i);
            end
        end
        [x, errorflag, errorcode] = solver(f, x, j1(j), j2(j), jacobian_flag, ...
                                           options.gstep, ...
                                           tolf, options.solve_tolx, maxit, ...
                                           options.trust_region_initial_step_bound_factor, ...
                                           options.debug, args{:});
        fre = true;
        if errorflag
            return
        end
    end
    fvec = feval(f, x, args{:});
    if max(abs(fvec))>tolf
        disp_verbose('Call solver on the full nonlinear problem.',options.verbosity)
        [x, errorflag, errorcode] = solver(f, x, 1:nn, 1:nn, jacobian_flag, ...
                                           options.gstep, tolf, options.solve_tolx, maxit, ...
                                           options.trust_region_initial_step_bound_factor, ...
                                           options.debug, args{:});
    end
    [fvec, fjac] = feval(f, x, args{:});
elseif options.solve_algo==3
    if jacobian_flag
        [x, errorcode] = csolve(f, x, f, tolf, maxit, args{:});
    else
        [x, errorcode] = csolve(f, x, [], tolf, maxit, args{:});
    end
    if errorcode==0
        errorflag = false;
    else
        errorflag = true;
    end
    [fvec, fjac] = feval(f, x, args{:});
elseif options.solve_algo==10
    % LMMCP
    olmmcp = options.lmmcp;
    [x, fvec, errorcode, ~, fjac] = lmmcp(f, x, olmmcp.lb, olmmcp.ub, olmmcp, args{:});
    if errorcode==1
        errorflag = false;
    else
        errorflag = true;
    end
    [fvec, fjac] = feval(f, x, args{:});
elseif options.solve_algo == 11
    % PATH mixed complementary problem
    % PATH linear mixed complementary problem
    if ~exist('mcppath')
        error(['PATH can''t be provided with Dynare. You need to install it ' ...
               'yourself and add its location to Matlab/Octave path before ' ...
               'running Dynare'])
    end
    omcppath = options.mcppath;
    global mcp_data
    mcp_data.func = f;
    mcp_data.args = args;
    try
        [x, fval, jac, mu] = pathmcp(x,omcppath.lb,omcppath.ub,'mcp_func',omcppath.A,omcppath.b,omcppath.t,omcppath.mu0);
    catch
        errorflag = true;
    end
    errorcode = nan; % There is no error code for this algorithm, as PATH is closed source it is unlikely we can fix that.
    [fvec, fjac] = feval(f, x, args{:});
elseif ismember(options.solve_algo, [13, 14])
    if ~jacobian_flag
        error('DYNARE_SOLVE: option solve_algo=13|14 needs computed Jacobian')
    end
    auxstruct = struct();
    if options.solve_algo == 14
        auxstruct.lhs = lhs;
        auxstruct.endo_names = endo_names;
        auxstruct.isloggedlhs = isloggedlhs;
        auxstruct.isauxdiffloggedrhs = isauxdiffloggedrhs;
    end
    [x, errorflag, errorcode] = block_trust_region(f, x, tolf, options.solve_tolx, maxit, options.trust_region_initial_step_bound_factor, options.debug, auxstruct, args{:});
    [fvec, fjac] = feval(f, x, args{:});
else
    error('DYNARE_SOLVE: option solve_algo must be one of [0,1,2,3,4,9,10,11,12,13,14]')
end
