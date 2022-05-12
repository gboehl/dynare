function [SE_values, Asympt_Var] = standard_errors(xparam, objective_function, Bounds, oo_, estim_params_, M_, options_mom_, Wopt_flag)
% [SE_values, Asympt_Var] = standard_errors(xparam, objective_function, Bounds, oo_, estim_params_, M_, options_mom_, Wopt_flag)
% -------------------------------------------------------------------------
% This function computes standard errors to the method of moments estimates
% Adapted from replication codes of
%  o Andreasen, Fernández-Villaverde, Rubio-Ramírez (2018): "The Pruned State-Space System for Non-Linear DSGE Models: Theory and Empirical Applications", Review of Economic Studies, 85(1):1-49.
% =========================================================================
% INPUTS
%   o xparam:                   value of estimated parameters as returned by set_prior()
%   o objective_function        string of objective function
%   o Bounds:                   structure containing parameter bounds
%   o oo_:                      structure for results
%   o estim_params_:            structure describing the estimated_parameters
%   o M_                        structure describing the model
%   o options_mom_:             structure information about all settings (specified by the user, preprocessor, and taken from global options_)
%   o Wopt_flag:                indicator whether the optimal weighting is actually used
% -------------------------------------------------------------------------
% OUTPUTS 
%   o SE_values                  [nparam x 1] vector of standard errors
%   o Asympt_Var                 [nparam x nparam] asymptotic covariance matrix
% -------------------------------------------------------------------------
% This function is called by
%  o mom.run.m
% -------------------------------------------------------------------------
% This function calls:
%  o get_the_name
%  o get_error_message
%  o mom.objective_function
%  o mom.optimal_weighting_matrix  
% =========================================================================
% Copyright (C) 2020-2021 Dynare Team
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
% -------------------------------------------------------------------------
% Author(s): 
% o Willi Mutschler (willi@mutschler.eu)
% o Johannes Pfeifer (jpfeifer@uni-koeln.de)
% =========================================================================

% Some dimensions
num_mom      = size(oo_.mom.model_moments,1);
dim_params   = size(xparam,1);
D            = zeros(num_mom,dim_params);
eps_value    = options_mom_.mom.se_tolx;

if strcmp(options_mom_.mom.mom_method,'GMM') && options_mom_.mom.analytic_standard_errors
    fprintf('\nComputing standard errors using analytical derivatives of moments\n');
    D = oo_.mom.model_moments_params_derivs; %already computed in objective function via get_perturbation_params.m
    idx_nan = find(any(isnan(D)));
    if any(idx_nan)
        for i = idx_nan            
             fprintf('No standard errors available for parameter %s\n',get_the_name(i,options_mom_.TeX, M_, estim_params_, options_mom_))
        end        
        warning('There are NaN in the analytical Jacobian of Moments. Check your bounds and/or priors, or use a different optimizer.')
        Asympt_Var = NaN(length(xparam),length(xparam));
        SE_values = NaN(length(xparam),1);
        return
    end
else    
    fprintf('\nComputing standard errors using numerical derivatives of moments\n');
    for i=1:dim_params
        %Positive step
        xparam_eps_p      = xparam;
        xparam_eps_p(i,1) = xparam_eps_p(i) + eps_value;
        [~, info_p, ~, ~,~, oo__p] = feval(objective_function, xparam_eps_p, Bounds, oo_, estim_params_, M_, options_mom_);

        % Negative step
        xparam_eps_m      = xparam;
        xparam_eps_m(i,1) = xparam_eps_m(i) - eps_value;
        [~, info_m,  ~, ~,~, oo__m] = feval(objective_function, xparam_eps_m, Bounds, oo_, estim_params_, M_, options_mom_);

        % The Jacobian:
        if nnz(info_p)==0 && nnz(info_m)==0
            D(:,i) = (oo__p.mom.model_moments - oo__m.mom.model_moments)/(2*eps_value);
        else
            problpar = get_the_name(i,options_mom_.TeX, M_, estim_params_, options_mom_);
            if info_p(1)==42
                warning('method_of_moments:info','Cannot compute the Jacobian using finite differences for parameter %s due to hitting the upper bound - no standard errors available.\n',problpar)
            else
                message_p = get_error_message(info_p, options_mom_);
            end
            if info_m(1)==41
                warning('method_of_moments:info','Cannot compute the Jacobian using finite differences for parameter %s due to hitting the lower bound - no standard errors available.\n',problpar)
            else
                message_m = get_error_message(info_m, options_mom_);        
            end
            if info_m(1)~=41 && info_p(1)~=42
                warning('method_of_moments:info','Cannot compute the Jacobian using finite differences for parameter %s - no standard errors available\n %s %s\nCheck your priors or use a different optimizer.\n',problpar, message_p, message_m)
            end
            Asympt_Var = NaN(length(xparam),length(xparam));
            SE_values = NaN(length(xparam),1);
            return
        end
    end
end

T = options_mom_.nobs; %Number of observations
if isfield(options_mom_,'variance_correction_factor')
    T = T*options_mom_.variance_correction_factor;
end

WW = oo_.mom.Sw'*oo_.mom.Sw;
if Wopt_flag
    % We have the optimal weighting matrix    
    Asympt_Var  = 1/T*((D'*WW*D)\eye(dim_params));
else
    % We do not have the optimal weighting matrix yet    
    WWopt      = mom.optimal_weighting_matrix(oo_.mom.m_data, oo_.mom.model_moments, options_mom_.mom.bartlett_kernel_lag);
    S          = WWopt\eye(size(WWopt,1));
    AA         = (D'*WW*D)\eye(dim_params);
    Asympt_Var = 1/T*AA*D'*WW*S*WW*D*AA;
end

SE_values   = sqrt(diag(Asympt_Var));
