function llik = likelihood(parameters)

% Evaluates the likelihood function.
%
% INPUTS
% - parameters    [char,double]    If row char array, possible values are 'posterior mode', 'posterior mean',
%                                  'posterior median', 'prior mode' or 'prior mean'. Otherwise, parmaters must
%                                  be a vector of doubles (arbitrary values for the parameters).
%
% OUTPUTS
% None

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

global M_ estim_params_ oo_ options_ bayestopt_

llik = evaluate_likelihood(parameters, M_, estim_params_, oo_, options_, bayestopt_);

if ~nargout
    dprintf('\nValue of the log-likelihood: %20.6f\n', llik);
    clear ('llik'); % Do not display the value returned by the function.
end