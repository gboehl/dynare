function [x, errorflag, errorcode] = newton_solve(func, x, jacobian_flag, gstep, tolf, tolx, maxit, solve_algo, varargin)

% Solves systems of non linear equations of several variables using a Newton solver, with three
% variants for the inner linear solver:
%  solve_algo = 6: use a sparse LU
%  solve_algo = 7: use GMRES
%  solve_algo = 8: use BiCGStab
%
% INPUTS
%    func:            name of the function to be solved
%    x:               guess values
%    jacobian_flag=true: jacobian given by the 'func' function
%    jacobian_flag=false: jacobian obtained numerically
%    gstep            increment multiplier in numerical derivative
%                     computation
%    tolf             tolerance for residuals
%    tolx             tolerance for solution variation
%    maxit            maximum number of iterations
%    solve_algo       from options_
%    varargin:        list of extra arguments to the function
%
% OUTPUTS
%    x:               results
%    errorflag=true:  the model can not be solved

% Copyright © 2024 Dynare Team
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

nn = length(x);

errorflag = false;

% Declaration needed for sharing with nested function
fvec = NaN(nn);
fjac = NaN(nn, nn);

function update_fvec_fjac
    if jacobian_flag
        [fvec, fjac] = feval(func, x, varargin{:});
    else
        fvec = feval(func, x, varargin{:});
        dh = max(abs(x),gstep(1)*ones(nn,1))*eps^(1/3);
        for j = 1:nn
            xdh = x;
            xdh(j) = xdh(j)+dh(j);
            t = feval(func, xdh, varargin{:});
            fjac(:, j) = (t - fvec) ./ dh(j);
        end
    end
end

update_fvec_fjac;

idInf = isinf(fvec);
idNan = isnan(fvec);
idCpx = ~isreal(fvec);

if any(idInf)
    disp('newton_solve: during the resolution of the non-linear system, the evaluation of the following equation(s) resulted in a non-finite number:')
    disp(idInf')
    errorcode = 0;
    errorflag = true;
    return
end

if any(idNan)
    disp('newton_solve: during the resolution of the non-linear system, the evaluation of the following equation(s) resulted in a nan:')
    disp(idNan')
    errorcode = 0;
    errorflag = true;
    return
end

if any(idNan)
    disp('newton_solve: during the resolution of the non-linear system, the evaluation of the following equation(s) resulted in a complex number:')
    disp(idCpx')
    errorcode = 0;
    errorflag = true;
    return
end

if max(abs(fvec)) < tolf
    % Initial guess is a solution
    errorcode = -1;
    return
end

for it = 1:maxit
    if solve_algo == 6
        p = fjac\fvec;
    else
        ilu_setup.type = 'ilutp';
        ilu_setup.droptol = 1e-10;
        [L1, U1] = ilu(fjac, ilu_setup);
        if solve_algo == 7
            p = gmres(fjac, fvec, [], [], [], L1, U1);
        elseif solve_algo == 8
            p = bicgstab(fjac, fvec, [], [], L1, U1);
        else
            error('newton_solve: invalid value for solve_algo')
        end
    end

    x = x - p;

    update_fvec_fjac;

    if max(abs(fvec)) < tolf
        errorcode = 1;
        return
    end
end

errorflag = true;
errorcode = 2;

end
