function oo_ = Jtest(xparam, objective_function, Woptflag, oo_, options_mom_, bayestopt_, BoundsInfo, estim_params_, M_, nobs)
% function oo_ = Jtest(xparam, objective_function, Woptflag, oo_, options_mom_, bayestopt_, BoundsInfo, estim_params_, M_, nobs)
% -------------------------------------------------------------------------
% Computes the J-test statistic and p-value for a GMM/SMM estimation
% =========================================================================
% INPUTS
%  xparam:              [vector]           estimated parameter vector
%  objective_function:  [function handle]  objective function
%  Woptflag:            [logical]          flag if optimal weighting matrix has already been computed
%  oo_:                 [struct]           results
%  options_mom_:        [struct]           options
%  bayestopt_:          [struct]           information on priors
%  BoundsInfo:          [struct]           bounds on parameters
%  estim_params_:       [struct]           information on estimated parameters
%  M_:                  [struct]           information on the model
%  nobs:                [scalar]           number of observations
% -------------------------------------------------------------------------
% OUTPUT
%  oo_:                 [struct]           updated results
% -------------------------------------------------------------------------
% This function is called by
%  o mom.run
% -------------------------------------------------------------------------
% This function calls
%  o mom.objective_function
%  o mom.optimal_weighting_matrix
% =========================================================================
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
% =========================================================================

if options_mom_.mom.mom_nbr > length(xparam)
    % Get optimal weighting matrix for J test, if necessary
    if ~Woptflag
        W_opt = mom.optimal_weighting_matrix(oo_.mom.m_data, oo_.mom.model_moments, options_mom_.mom.bartlett_kernel_lag);
        oo_J = oo_;
        oo_J.mom.Sw = chol(W_opt);
        fval = feval(objective_function, xparam, BoundsInfo, oo_J, estim_params_, M_, options_mom_);
    else
        fval = oo_.mom.Q;
    end
    % Compute J statistic
    if strcmp(options_mom_.mom.mom_method,'SMM')
        Variance_correction_factor = options_mom_.mom.variance_correction_factor;
    elseif strcmp(options_mom_.mom.mom_method,'GMM')
        Variance_correction_factor = 1;
    end
    oo_.mom.J_test.j_stat          = nobs*Variance_correction_factor*fval/options_mom_.mom.weighting_matrix_scaling_factor;
    oo_.mom.J_test.degrees_freedom = length(oo_.mom.model_moments)-length(xparam);
    oo_.mom.J_test.p_val           = 1-chi2cdf(oo_.mom.J_test.j_stat, oo_.mom.J_test.degrees_freedom);
    fprintf('\nValue of J-test statistic: %f\n',oo_.mom.J_test.j_stat);
    fprintf('p-value of J-test statistic: %f\n',oo_.mom.J_test.p_val);
end