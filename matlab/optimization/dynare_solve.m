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

% checking initial values
if jacobian_flag
    [fvec, fjac] = feval(f, x, varargin{:});
    wrong_initial_guess_flag = false;
    if ~all(isfinite(fvec)) || any(isinf(fjac(:))) || any(isnan((fjac(:)))) || any(~isreal(fvec)) || any(~isreal(fjac(:)))
        if ~ismember(options.solve_algo,[10,11]) && ~any(isnan(fvec)) && max(abs(fvec))< tolf
            % return if initial value solves the problem except if a mixed complementarity problem is to be solved (complementarity conditions may not be satisfied)
            % max([NaN, 0])=0, so explicitly exclude the case where fvec contains a NaN
            errorcode = -11;
            return;
        end
        if options.solve_randomize_initial_guess
            if any(~isreal(fvec)) || any(~isreal(fjac(:)))
                disp_verbose('dynare_solve: starting value results in complex values. Randomize initial guess...', options.verbosity)
            else
                disp_verbose('dynare_solve: starting value results in nonfinite/NaN value. Randomize initial guess...', options.verbosity)
            end
            % Let's try random numbers for the variables initialized with the default value.
            wrong_initial_guess_flag = true;
            % First try with positive numbers.
            tentative_number = 0;
            while wrong_initial_guess_flag && tentative_number<=in0*10
                tentative_number = tentative_number+1;
                x(idx) = rand(in0, 1)*10;
                [fvec, fjac] = feval(f, x, varargin{:});
                wrong_initial_guess_flag = ~all(isfinite(fvec)) || any(isinf(fjac(:))) || any(isnan((fjac(:)))) || any(~isreal(fvec)) || any(~isreal(fjac(:)));
            end
            % If all previous attempts failed, try with real numbers.
            tentative_number = 0;
            while wrong_initial_guess_flag && tentative_number<=in0*10
                tentative_number = tentative_number+1;
                x(idx) = randn(in0, 1)*10;
                [fvec, fjac] = feval(f, x, varargin{:});
                wrong_initial_guess_flag = ~all(isfinite(fvec)) || any(isinf(fjac(:))) || any(isnan((fjac(:)))) || any(~isreal(fvec)) || any(~isreal(fjac(:)));
            end
            % Last tentative, ff all previous attempts failed, try with negative numbers.
            tentative_number = 0;
            while wrong_initial_guess_flag && tentative_number<=in0*10
                tentative_number = tentative_number+1;
                x(idx) = -rand(in0, 1)*10;
                [fvec, fjac] = feval(f, x, varargin{:});
                wrong_initial_guess_flag = ~all(isfinite(fvec)) || any(isinf(fjac(:))) || any(isnan((fjac(:)))) || any(~isreal(fvec)) || any(~isreal(fjac(:)));
            end
        end
    end
else
    fvec = feval(f, x, varargin{:});
    fjac = zeros(nn, nn);
    if ~ismember(options.solve_algo,[10,11]) && ~any(isnan(fvec)) && max(abs(fvec)) < tolf
        % return if initial value solves the problem except if a mixed complementarity problem is to be solved (complementarity conditions may not be satisfied)
        % max([NaN, 0])=0, so explicitly exclude the case where fvec contains a NaN
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
            fvec = feval(f, x, varargin{:});
            wrong_initial_guess_flag = ~all(isfinite(fvec));
        end
        % If all previous attempts failed, try with real numbers.
        tentative_number = 0;
        while wrong_initial_guess_flag && tentative_number<=in0*10
            tentative_number = tentative_number+1;
            x(idx) = randn(in0, 1)*10;
            fvec = feval(f, x, varargin{:});
            wrong_initial_guess_flag = ~all(isfinite(fvec));
        end
        % Last tentative, ff all previous attempts failed, try with negative numbers.
        tentative_number = 0;
        while wrong_initial_guess_flag && tentative_number<=in0*10
            tentative_number = tentative_number+1;
            x(idx) = -rand(in0, 1)*10;
            fvec = feval(f, x, varargin{:});
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
    if isoctave
        options4fsolve = optimset('fsolve');
    else
        options4fsolve = optimoptions('fsolve');
    end
    if isoctave
        options4fsolve.MaxFunEvals = 50000;
        options4fsolve.MaxIter = maxit;
        options4fsolve.TolFun = tolf;
        options4fsolve.TolX = tolx;
        if jacobian_flag
            options4fsolve.Jacobian = 'on';
        else
            options4fsolve.Jacobian = 'off';
        end
    else
        options4fsolve.MaxFunctionEvaluations = 50000;
        options4fsolve.MaxIterations = maxit;
        options4fsolve.FunctionTolerance = tolf;
        options4fsolve.StepTolerance = tolx;
        options4fsolve.SpecifyObjectiveGradient = jacobian_flag;
    end
    %% NB: The Display option is accepted but not honoured under Octave (as of version 7)
    if options.debug
        options4fsolve.Display = 'final';
    else
        options4fsolve.Display = 'off';
    end
    %% This one comes last, so that the user can override Dynare
    if ~isempty(options.fsolve_options)
        if isoctave
            eval(['options4fsolve = optimset(options4fsolve,' options.fsolve_options ');']);
        else
            eval(['options4fsolve = optimoptions(options4fsolve,' options.fsolve_options ');']);
        end
    end

    if ~isoctave
        [x, fvec, errorcode, ~, fjac] = fsolve(f, x, options4fsolve, varargin{:});
    else
        % Under Octave, use a wrapper, since fsolve() does not have a 4th arg
        if ischar(f)
            f2 = str2func(f);
        else
            f2 = f;
        end
        [x, fvec, errorcode, ~, fjac] = fsolve(@(x) f2(x, varargin{:}), x, options4fsolve);
    end
    if errorcode==1
        errorflag = false;
    elseif errorcode>1
        if max(abs(fvec)) > tolf
            errorflag = true;
        else
            errorflag = false;
        end
    else
        errorflag = true;
    end
elseif ismember(options.solve_algo, [1, 12])
    %% NB: It is the responsibility of the caller to deal with the block decomposition if solve_algo=12
    [x, errorflag, errorcode] = solve1(f, x, 1:nn, 1:nn, jacobian_flag, options.gstep, tolf, tolx, maxit, [], options.debug, varargin{:});
    [fvec, fjac] = feval(f, x, varargin{:});
elseif options.solve_algo==9
    [x, errorflag, errorcode] = trust_region(f, x, 1:nn, 1:nn, jacobian_flag, options.gstep, tolf, tolx, maxit, options.trust_region_initial_step_bound_factor, options.debug, varargin{:});
    [fvec, fjac] = feval(f, x, varargin{:});
elseif ismember(options.solve_algo, [2, 4])
    if options.solve_algo == 2
        solver = @solve1;
    else
        solver = @trust_region;
    end
    if ~jacobian_flag
        fjac = zeros(nn,nn) ;
        dh = max(abs(x), options.gstep(1)*ones(nn,1))*eps^(1/3);
        for j = 1:nn
            xdh = x ;
            xdh(j) = xdh(j)+dh(j) ;
            fjac(:,j) = (feval(f, xdh, varargin{:})-fvec)./dh(j) ;
        end
    end
    [j1,j2,r,s] = dmperm(fjac);
    if options.debug
        disp(['DYNARE_SOLVE (solve_algo=2|4): number of blocks = ' num2str(length(r)-1)]);
    end
    for i=length(r)-1:-1:1
        blocklength = r(i+1)-r(i);
        j = r(i):r(i+1)-1;
        blockcolumns=s(i+1)-s(i);
        if blockcolumns ~= blocklength
            %non-square-block in DM; check whether initial value is solution
            [fval_check, fjac] = feval(f, x, varargin{:});
            if norm(fval_check(j1(j))) < tolf
                errorflag = false;
                errorcode = 0;
                continue
            end
        end
        if blockcolumns>=blocklength
            %(under-)determined block
            [x, errorflag, errorcode] = solver(f, x, j1(j), j2(j), jacobian_flag, ...
                options.gstep, ...
                tolf, options.solve_tolx, maxit, ...
                options.trust_region_initial_step_bound_factor, ...
                options.debug, varargin{:});
        else
            fprintf('\nDYNARE_SOLVE (solve_algo=2|4): the Dulmage-Mendelsohn decomposition returned a non-square block. This means that the Jacobian is singular. You may want to try another value for solve_algo.\n')
            %overdetermined block
            errorflag = true;
            errorcode = 0;
        end
        if errorflag
            return
        end
    end
    fvec = feval(f, x, varargin{:});
    if max(abs(fvec))>tolf
        disp_verbose('Call solver on the full nonlinear problem.',options.verbosity)
        [x, errorflag, errorcode] = solver(f, x, 1:nn, 1:nn, jacobian_flag, ...
                                           options.gstep, tolf, options.solve_tolx, maxit, ...
                                           options.trust_region_initial_step_bound_factor, ...
                                           options.debug, varargin{:});
    end
    [fvec, fjac] = feval(f, x, varargin{:});
elseif options.solve_algo==3
    if jacobian_flag
        [x, errorcode] = csolve(f, x, f, tolf, maxit, varargin{:});
    else
        [x, errorcode] = csolve(f, x, [], tolf, maxit, varargin{:});
    end
    if errorcode==0
        errorflag = false;
    else
        errorflag = true;
    end
    [fvec, fjac] = feval(f, x, varargin{:});
elseif options.solve_algo==10
    % LMMCP
    olmmcp = options.lmmcp;
    [x, fvec, errorcode, ~, fjac] = lmmcp(f, x, olmmcp.lb, olmmcp.ub, olmmcp, varargin{:});
    eq_to_check=find(isfinite(olmmcp.lb) | isfinite(olmmcp.ub));
    eq_to_ignore=eq_to_check(x(eq_to_check,:)<=olmmcp.lb(eq_to_check)+eps | x(eq_to_check,:)>=olmmcp.ub(eq_to_check)-eps);
    fvec(eq_to_ignore)=0;
    if errorcode==1
        errorflag = false;
    else
        errorflag = true;
    end
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
    mcp_data.args = varargin;
    try
        x = pathmcp(x,omcppath.lb,omcppath.ub,'mcp_func',omcppath.A,omcppath.b,omcppath.t,omcppath.mu0);
    catch
        errorflag = true;
    end
    errorcode = nan; % There is no error code for this algorithm, as PATH is closed source it is unlikely we can fix that.
    eq_to_check=find(isfinite(omcppath.lb) | isfinite(omcppath.ub));
    eq_to_ignore=eq_to_check(x(eq_to_check,:)<=omcppath.lb(eq_to_check)+eps | x(eq_to_check,:)>=omcppath.ub(eq_to_check)-eps);
    fvec(eq_to_ignore)=0;
elseif ismember(options.solve_algo, [13, 14])
    %% NB: It is the responsibility of the caller to deal with the block decomposition if solve_algo=14
    if ~jacobian_flag
        error('DYNARE_SOLVE: option solve_algo=13 needs computed Jacobian')
    end
    [x, errorflag, errorcode] = block_trust_region(f, x, tolf, options.solve_tolx, maxit, ...
                                                   options.trust_region_initial_step_bound_factor, ...
                                                   options.solve_algo == 13, ... % Only block-decompose with Dulmage-Mendelsohn for 13, not for 14
                                                   options.debug, varargin{:});
    [fvec, fjac] = feval(f, x, varargin{:});
else
    error('DYNARE_SOLVE: option solve_algo must be one of [0,1,2,3,4,9,10,11,12,13,14]')
end
