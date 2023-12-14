function [llik,parameters] = evaluate_likelihood(parameters,M_,estim_params_,oo_,options_,bayestopt_)
% [llik,parameters] = evaluate_likelihood(parameters,M_,estim_params_,oo_,options_,bayestopt_)
% Evaluate the logged likelihood at parameters.
%
% INPUTS
%    o parameters  a string ('posterior mode','posterior mean','posterior median','prior mode','prior mean') or a vector of values for
%                  the (estimated) parameters of the model.
%    o M_               [structure]  Definition of the model
%    o estim_params_    [structure] characterizing parameters to be estimated
%    o oo_              [structure]  Storage of results
%    o options_         [structure]  Options
%    o bayestopt_       [structure]  describing the priors
%
% OUTPUTS
%    o ldens            [double]  value of the sample logged density at parameters.
%    o parameters       [double]  vector of values for the estimated parameters.
%
% SPECIAL REQUIREMENTS
%    None
%
% REMARKS
% [1] This function cannot evaluate the likelihood of a dsge-var model...
% [2] This function use persistent variables for the dataset_ and the description of the missing observations. Consequently, if this function
%     is called more than once (by changing the value of parameters) the sample *must not* change.

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

persistent dataset_ dataset_info

if nargin==0
    parameters = 'posterior mode';
end

if ischar(parameters)
    switch parameters
        case 'posterior mode'
            parameters = get_posterior_parameters('mode',M_,estim_params_,oo_,options_);
        case 'posterior mean'
            parameters = get_posterior_parameters('mean',M_,estim_params_,oo_,options_);
        case 'posterior median'
            parameters = get_posterior_parameters('median',M_,estim_params_,oo_,options_);
        case 'prior mode'
            parameters = bayestopt_.p5(:);
        case 'prior mean'
            parameters = bayestopt_.p1;
        otherwise
            disp('eval_likelihood:: If the input argument is a string, then it has to be equal to:')
            disp('                   ''posterior mode'', ')
            disp('                   ''posterior mean'', ')
            disp('                   ''posterior median'', ')
            disp('                   ''prior mode'' or')
            disp('                   ''prior mean''.')
            error
    end
end

if isempty(dataset_)
    [dataset_, dataset_info] = makedataset(options_);
end
options_=select_qz_criterium_value(options_);

if ~isempty(bayestopt_) && any(bayestopt_.pshape > 0)
    % Plot prior densities.
    % Set prior bounds
    bounds = prior_bounds(bayestopt_, options_.prior_trunc);
else  % estimated parameters but no declared priors
    % No priors are declared so Dynare will estimate the model by
    % maximum likelihood with inequality constraints for the parameters.
    [~,~,~,lb,ub] = set_prior(estim_params_,M_,options_);
    bounds.lb = lb;
    bounds.ub = ub;
end


if options_.occbin.likelihood.status && options_.occbin.likelihood.inversion_filter
    llik = -occbin.IVF_posterior(parameters,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_.dr, oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
else
    llik = -dsge_likelihood(parameters,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_.dr, oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
end
ldens = evaluate_prior(parameters,M_,estim_params_,oo_,options_,bayestopt_);
llik = llik - ldens;