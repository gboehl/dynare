function [tlogpostkernel,loglikelihood] = tempered_likelihood(postkernelfun, xparam, lambda, Prior)

% Evaluate tempered likelihood (posterior kernel)
%
% INPUTS
% - postkernelfun       [handle]   Function handle for the opposite of the  posterior kernel.
% - xparam              [double]   n×1 vector of parameters.
% - lambda              [double]   scalar between 0 and 1, weight on the posterior kernel.
% - Prior               [dprior]   Prior specification.
%
% OUTPUTS
% - tlogpostkernel      [double]   scalar, value of the tempered posterior kernel.
% - loglikelihood       [double]   scalar, value of the log likelihood.

% Copyright © 2022-2023 Dynare Team
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

logpostkernel = -postkernelfun(xparam);
logprior = Prior.density(xparam);
loglikelihood = logpostkernel-logprior;
tlogpostkernel = lambda*loglikelihood + logprior;
