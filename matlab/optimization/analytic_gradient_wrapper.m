function [fval, grad, hess, exit_flag]=analytic_gradient_wrapper(x, fcn, varargin)
%function [fval, grad, hess, exitflag]=analytic_gradient_wrapper(x, fcn, varargin)
% Encapsulates an objective function to be minimized for use with Matlab
% optimizers
%
% INPUTS
% - x             [double]    n*1 vector of instrument values.
% - fcn           [fhandle]   objective function.
% - varagin       [cell]      additional parameters for fcn.
%
% OUTPUTS
% - fval          [double]    scalar, value of the objective function at x.
% - grad                      gradient of the objective function
% - hess                      Hessian of the objective function
% - exit_flag     [integer]   scalar, flag returned by

% Copyright © 2021 Dynare Team
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

[fval, ~, exit_flag, grad, hess] = fcn(x, varargin{:});
if size(grad,2)==1
    grad=grad'; %should be row vector for Matlab; exception lsqnonlin where Jacobian is required
end