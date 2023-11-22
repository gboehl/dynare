function [x, errorflag, info] = trust_region(objfun, x, j1, j2, jacobianflag, gstep, tolf, tolx, maxiter, factor, debug, varargin)

% Solves systems of non linear equations of several variables, using a
% trust-region method.
%
% INPUTS
% - objfun         [function handle, char]      name of the routine evaluating the system of nonlinear equations (and possibly jacobian).
% - x              [double]                     n×1 vector, initial guess for the solution.
% - j1             [integer]                    vector, equation indices defining a subproblem to be solved.
% - j2             [integer]                    vector, unknown variable indices to be solved for in the subproblem.
% - jacobianflag   [logical]                    scalar, if true the jacobian matrix is expected to be returned as a second output argument when calling objfun, otherwise
%                                               the jacobian is computed numerically.
% - gstep          [double]                     scalar, increment multiplier in numerical derivative computation (only used if jacobianflag value is false).
% - tolf           [double]                     scalar, tolerance for residuals.
% - tolx           [double]                     scalar, tolerance for solution variation.
% - maxiter        [integer]                    scalar, maximum number of iterations;
% - factor         [double]                     scalar, determines the initial step bound.
% - debug          [logical]                    scalar, dummy argument.
% - varargin:      [cell]                       list of additional arguments to be passed to objfun.
%
% OUTPUTS
% - x              [double]                     n⨱1 vector, solution of the nonlinear system of equations.
% - errorflag      [logical]                    scalar, false iff nonlinear solver is successful.
% - info           [integer]                    scalar, information about the failure.
%
% REMARKS
% [1] j1 and j2 muyst have the same number of elements.
% [2] debug is here for compatibility purpose (see solve1), it does not affect the output.
% [3] Possible values for info are:
%
%       -1 if the initial guess is a solution of the nonlinear system of equations.
%        0 if the nonlinear solver failed because the problem is ill behaved at the initial guess.
%        1 if the nonlinear solver found a solution.
%        2 if the maximum number of iterations has been reached.
%        3 if spurious convergence (trust region radius is too small).
%        4 if iteration is not making good progress, as measured by the improvement from the last 15 iterations.
%        5 if no further improvement in the approximate solution x is possible (xtol is too small).

% Copyright © 2014-2023 Dynare Team
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

% Convert to function handle if necessary
if ischar(objfun)
    objfun = str2func(objfun);
end

%
% Set constants
%

radiusfactor = .5;
actualreductionthreshold = .001;
relativereductionthreshold = .1;

% number of equations
n = length(j1);

%
% Initialization
%

errorflag = true;

iter = 1;
info = 0;
delta = 0.0;

ncsucc = 0; ncslow = 0;

xinit = x;

%
% Attempt to evaluate the residuals and jacobian matrix on the initial guess
%

try
    if jacobianflag
        [fval, fjac] = objfun(x, varargin{:});
        fval = fval(j1);
        fjac = fjac(j1,j2);
    else
        fval = objfun(x, varargin{:});
        fval = fval(j1);
        dh = max(abs(x(j2)), gstep(1)*ones(n,1))*eps^(1/3);
        fjac = zeros(n);
        for j = 1:n
            xdh = x;
            xdh(j2(j)) = xdh(j2(j))+dh(j);
            Fval = objfun(xdh, varargin{:});
            fjac(:,j) = (Fval(j1)-fval)./dh(j);
        end
    end
    fnorm = norm(fval);
catch
    % System of equation cannot be evaluated at the initial guess.
    return
end

if any(isnan(fval)) || any(isinf(fval)) || any(~isreal(fval)) || any(isnan(fjac(:))) || any(isinf(fjac(:))) || any(~isreal(fjac(:)))
    % System of equations is ill-behaved at the initial guess.
    return
end

if fnorm<tolf
    % Initial guess is a solution.
    errorflag = false;
    info = -1;
    return
end

%
% Main loop
%

while iter<=maxiter && ~info
    % Compute the columns norm ofr the Jacobian matrix.
    fjacnorm = transpose(sqrt(sum((fjac.*fjac))));
    if iter==1
        % On the first iteration, calculate the norm of the scaled vector of unknowns x
        % and initialize the step bound delta. Scaling is done according to the norms of
        % the columns of the initial jacobian.
        fjacnorm__ = fjacnorm;
        fjacnorm__(fjacnorm<eps(1.0)) = 1.0;
        xnorm = norm(fjacnorm__.*x(j2));
        if xnorm>0
            delta = xnorm*factor;
        else
            delta = factor;
        end
    else
        fjacnorm__ = max(.1*fjacnorm__, fjacnorm);
        xnorm = norm(fjacnorm__.*x(j2));
    end
    % Determine the direction p (with trust region model defined in dogleg routine).
    p = dogleg(fjac, fval, fjacnorm__, delta);
    % Compute the norm of p.
    pnorm = norm(fjacnorm__.*p);
    xx = x;
    x0 = x;
    % Move along the direction p. Set a candidate value for x and predicted improvement for f.
    xx(j2) = x(j2) - p;
    ww = fval - fjac*p;
    % Evaluate the function at xx and calculate its norm.
    try
        fval1 = objfun(xx, varargin{:});
        fval1 = fval1(j1);
    catch
        % If evaluation of the residuals returns an error, then restart but with a smaller radius of the trust region.
        delta = delta*radiusfactor;
        iter = iter+1;
        continue
    end
    if any(isnan(fval1)) || any(isinf(fval1)) || any(~isreal(fval1))
        % If evaluation of the residuals returns a NaN, an infinite number or a complex number, then restart but with a smaller radius of the trust region.
        delta = delta*radiusfactor;
        iter = iter+1;
        continue
    end
    fnorm1 = norm(fval1);
    if fnorm1<tolf
        x = xx;
        errorflag = false;
        info = 1;
        continue
    end
    % Compute the scaled actual reduction.
    if fnorm1<fnorm
        actualreduction = 1.0-(fnorm1/fnorm)^2;
    else
        actualreduction = -1.0;
    end
    % Compute the scaled predicted reduction and the ratio of the actual to the
    % predicted reduction.
    tt = norm(ww);
    if tt<fnorm
        predictedreduction = 1.0 - (tt/fnorm)^2;
    else
        predictedreduction = 0.0;
    end
    if predictedreduction>0
        ratio = actualreduction/predictedreduction;
    else
        ratio = 0;
    end
    % Update the radius of the trust region if need be.
    if iter==1
        % On first iteration adjust the initial step bound.
        delta = min(delta, pnorm);
    end
    if ratio<relativereductionthreshold
        % Reduction is much smaller than predicted… Reduce the radius of the trust region.
        ncsucc = 0;
        delta = delta*radiusfactor;
    else
        ncsucc = ncsucc + 1;
        if abs(ratio-1.0)<relativereductionthreshold
            delta = pnorm/radiusfactor;
        elseif ratio>=radiusfactor || ncsucc>1
            delta = max(delta, pnorm/radiusfactor);
        end
    end
    if ratio>1e-4
        % Successful iteration. Update x, xnorm, fval, fnorm and fjac.
        x = xx;
        fval = fval1;
        xnorm = norm(fjacnorm__.*x(j2));
        fnorm = fnorm1;
    end
    % Determine the progress of the iteration.
    ncslow = ncslow+1;
    if actualreduction>=actualreductionthreshold
        ncslow = 0;
    end
    iter = iter+1;
    if iter==maxiter
        info = 2;
        x = xinit;
        continue
    end
    if delta<tolx*xnorm
        info = 3;
        x = xinit;
        errorflag = true;
        continue
    end
    % Tests for termination and stringent tolerances.
    if max(.1*delta, pnorm)<=10*eps(xnorm)*xnorm
        % xtol is too small. no further improvement in
        % the approximate solution x is possible.
        info = 5;
        x = xinit;
        errorflag = true;
        continue
    end
    if ncslow==15
        info = 4;
        x = xinit;
        errorflag = true;
        continue
    end
    % Compute the jacobian for the next iteration.
    fjac0 = fjac;
    if jacobianflag
        try
            [~, fjac] = objfun(x, varargin{:});
            fjac = fjac(j1,j2);
        catch
            % If evaluation of the Jacobian matrix returns an error, then restart but with a smaller radius of the trust region.
            x = x0;
            fjac = fjac0;
            delta = delta*radiusfactor;
            continue
        end
    else
        dh = max(abs(x(j2)), gstep(1)*ones(n,1))*eps^(1/3);
        for j = 1:n
            xdh = x;
            xdh(j2(j)) = xdh(j2(j))+dh(j);
            Fval = objfun(xdh, varargin{:});
            fjac(:,j) = (Fval(j1)-fval)./dh(j);
        end
    end
    if any(isnan(fjac(:))) || any(isinf(fjac(:))) || any(~isreal(fjac(:)))
        % If evaluation of the Jacobian matrix returns NaNs, an infinite numbers or a complex numbers, then restart but with a smaller radius of the trust region.
        x = x0;
        fjac = fjac0;
        delta = delta*radiusfactor;
    end
end


function x = dogleg (r, b, d, delta)
% Compute the Gauss-Newton direction.
if isoctave
    % The decomposition() function does not exist in Octave
    x = r \ b;
else
    x = decomposition(r, 'CheckCondition', false) \ b;
end
% Compute norm of scaled x
qnorm = norm(d.*x);
if qnorm<=delta
    % Gauss-Newton direction is acceptable. There is nothing to do here.
else
    % Gauss-Newton direction is not acceptable…
    % Compute the scale gradient direction and its norm
    s = (r'*b)./d;
    gnorm = norm(s);
    if gnorm>0
        % Normalize and rescale → gradient direction.
        s = (s/gnorm)./d;
        % Get the line minimizer in s direction.
        temp0 = norm(r*s);
        sgnorm = gnorm/(temp0*temp0);
        if sgnorm<delta
            % The scaled gradient direction is not acceptable…
            % Compute the point along the dogleg at which the
            % quadratic is minimized.
            bnorm = norm(b);
            temp1 = delta/qnorm;
            temp2 = sgnorm/delta;
            temp0 = bnorm*bnorm*temp2/(gnorm*qnorm);
            temp0 = temp0 - temp1*temp2*temp2 + sqrt((temp0-temp1)^2+(1.0-temp1*temp1)*(1.0-temp2*temp2));
            alpha = temp1*(1.0-temp2*temp2)/temp0;
        else
            % The scaled gradient direction is acceptable.
            alpha = 0.0;
        end
    else
        % If the norm of the scaled gradient direction is zero.
        alpha = delta/qnorm;
        sgnorm = 0.0;
    end
    % Form the appropriate  convex combination of the Gauss-Newton direction and the
    % scaled gradient direction.
    if alpha>0
        x = alpha*x + (1.0-alpha)*min(sgnorm, delta)*s;
    else %prevent zero weight on Inf evaluating to NaN
        x = min(sgnorm, delta)*s;
    end
end
