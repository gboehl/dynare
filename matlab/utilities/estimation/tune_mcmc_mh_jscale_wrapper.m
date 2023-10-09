function mh_jscale = tune_mcmc_mh_jscale_wrapper(invhess, options_, M_, objective_function, xparam1, bounds, varargin)
% function mh_jscale = tune_mcmc_mh_jscale_wrapper(invhess, options_, M_, objective_function, xparam1, bounds, varargin)
% -------------------------------------------------------------------------
% Wrapper to call the algorithm to tune the jumping scale parameter for the
% Metropolis-Hastings algorithm; currently only supports RW-MH algorithm.
% -------------------------------------------------------------------------
% INPUTS
%  o invhess:                 [matrix] jumping covariance matrix
%  o options_:                [struct] Dynare options
%  o M_:                      [struct] Dynare model structure
%  o objective_function:      [function handle] objective function
%  o xparam1:                 [vector] vector of estimated parameters at the mode
%  o bounds:                  [struct] structure containing information on bounds
%  o varargin:                [cell] additional arguments to be passed to the objective function
% -------------------------------------------------------------------------
% OUTPUTS
%  o mh_jscale:               [double] tuned jumping scale parameter
% -------------------------------------------------------------------------
% This function is called by
%  o dynare_estimation_1
%  o mom.run
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

posterior_sampler_options_temp = options_.posterior_sampler_options.current_options;
posterior_sampler_options_temp.invhess = invhess;
posterior_sampler_options_temp = check_posterior_sampler_options(posterior_sampler_options_temp, M_.fname, M_.dname, options_);
opt = options_.mh_tune_jscale;
opt.rwmh = options_.posterior_sampler_options.rwmh;
mh_jscale = calibrate_mh_scale_parameter(objective_function, ...
                                          posterior_sampler_options_temp.invhess, xparam1, [bounds.lb,bounds.ub], ...
                                          opt, varargin{:});