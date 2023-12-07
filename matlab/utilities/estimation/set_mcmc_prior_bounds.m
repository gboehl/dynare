function bounds = set_mcmc_prior_bounds(xparam, bayestopt_, options_, stringForErrors)
% function bounds = set_mcmc_prior_bounds(xparam, bayestopt_, options_, stringForErrors)
% -------------------------------------------------------------------------
% Reset bounds as lower and upper bounds must only be operational during mode-finding
% =========================================================================
% INPUTS
%  o xparam:                 [vector] vector of parameters
%  o bayestopt_:             [struct] information on priors
%  o options_:               [struct] Dynare options
%  o stringForErrors:        [string] string to be used in error messages
% -------------------------------------------------------------------------
% OUTPUTS
%  o bounds:                 [struct] structure with fields lb and ub
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

bounds = prior_bounds(bayestopt_, options_.prior_trunc);
outside_bound_pars = find(xparam < bounds.lb | xparam > bounds.ub);
if ~isempty(outside_bound_pars)
    for ii = 1:length(outside_bound_pars)
        outside_bound_par_names{ii,1} = get_the_name(ii,0,M_,estim_params_,options_.varobs);
    end
    disp_string = [outside_bound_par_names{1,:}];
    for ii = 2:size(outside_bound_par_names,1)
        disp_string = [disp_string,', ',outside_bound_par_names{ii,:}];
    end
    if options_.prior_trunc > 0
        error(['%s: Mode value(s) of ', disp_string ,' are outside parameter bounds. Potentially, you should set prior_trunc=0!'],stringForErrors);
    else
        error(['%s: Mode value(s) of ', disp_string ,' are outside parameter bounds!'],stringForErrors);
    end
end