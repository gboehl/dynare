function [xparams, lpd, hessian_mat] = ...
    maximize_prior_density(iparams, names, options_, M_, Prior, estim_params_, oo_)

% Maximizes the logged prior density using Chris Sims' optimization routine.
%
% INPUTS
% - iparams                        [double]   vector of initial parameters.
% - Prior              [dprior]   vector specifying prior densities shapes.
% - DynareOptions      [struct]   Options, AKA options_
% - DynareModel        [struct]   Model description, AKA M_
% - EstimatedParams    [struct]   Info about estimated parameters, AKA estimated_params_
% - DynareResults      [struct]   Results, AKA oo_

%
%
%   prior_shape                    [integer]  vector specifying prior densities shapes.
%   prior_hyperparameter_1         [double]   vector, first hyperparameter.
%   prior_hyperparameter_2         [double]   vector, second hyperparameter.
%   prior_inf_bound                [double]   vector, prior's lower bound.
%   prior_sup_bound                [double]   vector, prior's upper bound.
%   options_                       [structure] describing the options
%   bayestopt_                     [structure] describing the priors
%   M_                             [structure] describing the model
%   estim_params_                  [structure] characterizing parameters to be estimated
%   oo_                            [structure] storing the results
%
% OUTPUTS
%   xparams       [double]  vector, prior mode.
%   lpd           [double]  scalar, value of the logged prior density at the mode.
%   hessian_mat   [double]  matrix, Hessian matrix at the prior mode.

% Copyright © 2009-2023 Dynare Team
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

[xparams, lpd, ~, hessian_mat] = dynare_minimize_objective('minus_logged_prior_density', ...
                                                                  iparams, ...
                                                                  options_.mode_compute, ...
                                                                  options_, ...
                                                                  [Prior.p3, Prior.p4], ...
                                                                  names, ...
                                                                  [], ...
                                                                  [], ...
                                                                  Prior, ...
                                                                  options_, ...
                                                                  M_, ...
                                                                  estim_params_, ...
                                                                  oo_);

lpd = -lpd;
