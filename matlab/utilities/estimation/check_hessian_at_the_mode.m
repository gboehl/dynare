function check_hessian_at_the_mode(hessian_xparam1, xparam1, M_, estim_params_, options_, bounds)
% check_hessian_at_the_mode(hessian_xparam1, xparam1, M_, estim_params_, options_, bounds)
% -------------------------------------------------------------------------
% This function checks whether the hessian matrix at the mode is positive definite.
% -------------------------------------------------------------------------
% INPUTS
%  o hessian_xparam1:        [matrix] hessian matrix at the mode
%  o xparam1:                [vector] vector of parameter values at the mode
%  o M_:                     [structure] information about model
%  o estim_params_:          [structure] information about estimated parameters
%  o options_:               [structure] information about options
%  o bounds:                 [structure] information about bounds
% -------------------------------------------------------------------------
% OUTPUTS
%   none, displays a warning message if the hessian matrix is not positive definite
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

try
    chol(hessian_xparam1);
catch
    tol_bounds = 1.e-10;
    skipline()
    disp('OPTIMIZATION PROBLEM!')
    disp(' (minus) the hessian matrix at the "mode" is not positive definite!')
    disp('=> variance of the estimated parameters are not positive.')
    disp('You should try to change the initial values of the parameters using')
    disp('the estimated_params_init block, or use another optimization routine.')
    params_at_bound = find(abs(xparam1-bounds.ub)<tol_bounds | abs(xparam1-bounds.lb)<tol_bounds);
    if ~isempty(params_at_bound)
        for ii=1:length(params_at_bound)
            params_at_bound_name{ii,1}=get_the_name(params_at_bound(ii),0,M_,estim_params_,options_.varobs);
        end
        disp_string=[params_at_bound_name{1,:}];
        for ii=2:size(params_at_bound_name,1)
            disp_string=[disp_string,', ',params_at_bound_name{ii,:}];
        end
        fprintf('\nThe following parameters are at the bound: %s\n', disp_string)
        fprintf('Some potential solutions are:\n')
        fprintf('   - Check your model for mistakes.\n')
        fprintf('   - Check whether model and data are consistent (correct observation equation).\n')
        fprintf('   - Shut off prior_trunc.\n')
        fprintf('   - Change the optimization bounds.\n')
        fprintf('   - Use a different mode_compute like 6 or 9.\n')
        fprintf('   - Check whether the parameters estimated are identified.\n')
        fprintf('   - Check prior shape (e.g. Inf density at bound(s)).\n')
        fprintf('   - Increase the informativeness of the prior.\n')
    end
    warning('The results below are most likely wrong!');
end