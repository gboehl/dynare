function [xparams, logpost] = GetOneDraw(distribution, M_, estim_params_, oo_, options_, bayestopt_)

% draws one parameter vector and its posterior from MCMC or the prior
%
% INPUTS
%    type:      [string]       'posterior': draw from MCMC draws
%                              'prior': draw from prior
%    M_         [structure]     Definition of the model
%    estim_params_ [structure]  characterizing parameters to be estimated
%    oo_         [structure]    Storage of results
%    options_    [structure]    Options
%    bayestopt_  [structure]    describing the priors
%
% OUTPUTS
%    xparams:   vector of estimated parameters (drawn from posterior or prior distribution)
%    logpost:   log of the posterior density of this parameter vector
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2005-2023 Dynare Team
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

if isprior(distribution)
    xparams = distribution.draw();
    if nargout>1
        logpost = evaluate_posterior_kernel(xparams, M_, estim_params_, oo_, options_, bayestopt_);
    end
elseif ischar(distribution) && strcmpi(distribution, 'posterior')
    [xparams, logpost] = metropolis_draw(0);
else
    error('GetOneDraw:: Wrong inputs.')
end
