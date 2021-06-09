function [x, objectivevalue, errorflag] = gauss_newton(fun, x0)

% Minimization of the sum of squared residuals with the Gauss-Newton algorithm.
%
% The objective is to minimize:
%
%        fun(x)'*fun(x)
%
% with respect to x.
%
% INPUTS:
% - funres          [handle]   Function from Rᵖ to Rⁿ, which given parameters (x) return the residuals of a non linear equation.
% - x0              [double]   1×p vector, initial guess.
%
% OUTPUTS:
% - x               [double]   1×p vector, vector of parameters minimizing the sum of squared residuals.
% - objectivevalue  [double]   scalar, optimal value of the objective.
% - errorflag       [integer]  scalar, nonzero if algorithm did not converge.

% Copyright (C) 2018 Dynare Team
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

maxIter = 100;           % Maximum number of iteration.
h = sqrt(eps(1.0));      % Pertubation size for numerical computation of the Jacobian matrix
xtol = 1e-6;             % Stopping criterion.
phi = .5*(sqrt(5)+1.0);  % Golden number.

errorflag = 0;

noconvergence = true;
counter = 0;

while noconvergence
    % Compute residuals and descent direction
    [r0, J] = jacobian(fun, x0, h);
    d = pinv(J)*r0;
    % Update parameters
    x1 = x0+d;
    % Test if the step actually reduce the sum of squared residuals.
    r1 = fun(x1);
    s0 = r0'*r0;
    s1 = r1'*r1;
    if s1>s0
        % Gauss-Newton step increased the Sum of Squared Residuals...
        % We search for another point in the same direction using Golden section search.
        l1 = 0;
        l2 = 1;
        L1 = l2-(l2-l1)/phi;
        L2 = l1+(l2-l1)/phi;
        while abs(L1-L2)>1e-6
            if ssr(x0+L1*d)<ssr(x0+L2*d)
                l2 = L2;
            else
                l1 = L1;
            end
            L1 = l2-(l2-l1)/phi;
            L2 = l1+(l2-l1)/phi;
        end
        scale = .5*(l1+l2);
        x1 = x0+scale*d;
    else
        scale = 1.0;
    end
    noconvergence = max(abs(x1-x0))>xtol;
    counter = counter+1;
    x0 = x1;
    if counter>maxIter
        break
    end
end

x = x0;
objectivevalue = s1;

errorflag = isequal(counter, maxIter+1);

function s = ssr(x)
    % Evaluate the sum of square residuals.
    r = fun(x);
    s = r'*r;
end

end