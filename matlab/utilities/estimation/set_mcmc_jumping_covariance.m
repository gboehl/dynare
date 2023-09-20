function invhess = set_mcmc_jumping_covariance(invhess, xparam_nbr, MCMC_jumping_covariance, bayestopt_, stringForErrors)
% function invhess = set_mcmc_jumping_covariance(invhess, xparam_nbr, MCMC_jumping_covariance, bayestopt_, stringForErrors)
% -------------------------------------------------------------------------
% sets the jumping covariance matrix for the MCMC algorithm
% -------------------------------------------------------------------------
% INPUTS
%  o invhess:                 [matrix] already computed inverse of the hessian matrix
%  o xparam_nbr:              [integer] number of estimated parameters
%  o MCMC_jumping_covariance: [string] name of option or file setting the jumping covariance matrix
%  o bayestopt_:              [struct] information on priors
%  o stringForErrors:         [string] string to be used in error messages
% -------------------------------------------------------------------------
% OUTPUTS
%  o invhess:                 [matrix] jumping covariance matrix
% -------------------------------------------------------------------------
% This function is called by
%  o dynare_estimation_1.m
% -------------------------------------------------------------------------

% Copyright © 2023 Dynare Team
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

switch MCMC_jumping_covariance
    case 'hessian' % do nothing and use hessian from previous mode optimization
    case 'prior_variance' % use prior variance
        if any(isinf(bayestopt_.p2))
            error('%s: Infinite prior variances detected. You cannot use the prior variances as the proposal density, if some variances are Inf.',stringForErrors);
        else
            hess = diag(1./(bayestopt_.p2.^2));
        end
        hsd = sqrt(diag(hess));
        invhess = inv(hess./(hsd*hsd'))./(hsd*hsd');
    case 'identity_matrix' % use identity
        invhess = eye(xparam_nbr);
    otherwise % user specified matrix in file
        try
            load(MCMC_jumping_covariance,'jumping_covariance')
            hess = jumping_covariance;
        catch
            error(['%s: No matrix named ''jumping_covariance'' could be found in ',options_.MCMC_jumping_covariance,'.mat!'],stringForErrors);
        end
        [nrow, ncol] = size(hess);
        if ~isequal(nrow,ncol) && ~isequal(nrow,xparam_nbr) % check if square and right size
            error(['%s: ''jumping_covariance'' matrix (loaded from ',options_.MCMC_jumping_covariance,'.mat) must be square and have ',num2str(xparam_nbr),' rows and columns!'],stringForErrors);
        end
        try % check for positive definiteness
            chol(hess);
            hsd = sqrt(diag(hess));
            invhess = inv(hess./(hsd*hsd'))./(hsd*hsd');
        catch
            error('%s: Specified ''MCMC_jumping_covariance'' is not positive definite!',stringForErrors);
        end
end