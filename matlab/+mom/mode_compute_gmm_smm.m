function [xparam1, oo_, Woptflag] = mode_compute_gmm_smm(xparam0, objective_function, oo_, M_, options_mom_, estim_params_, bayestopt_, BoundsInfo)
% function [xparam1, oo_, Woptflag] = mode_compute_gmm_smm(xparam0, objective_function, oo_, M_, options_mom_, estim_params_, bayestopt_, BoundsInfo)
% -------------------------------------------------------------------------
% Iterated method of moments for GMM and SMM, computes the minimum of the
% objective function (distance between data moments and model moments)
% for a sequence of optimizers and GMM/SMM iterations with different
% weighting matrices.
% =========================================================================
% INPUTS
% xparam0:             [vector]       vector of initialized parameters
% objective_function:  [func handle]  name of the objective function
% oo_:                 [structure]    results
% M_:                  [structure]    model information
% options_mom_:        [structure]    options
% estim_params_:       [structure]    information on estimated parameters
% bayestopt_:          [structure]    information on priors
% BoundsInfo:          [structure]    bounds for optimization
% -------------------------------------------------------------------------
% OUTPUT
% xparam1:             [vector]       mode of objective function
% oo_:                 [structure]    updated results
% Woptflag:            [logical]      true if optimal weighting matrix was computed
% -------------------------------------------------------------------------
% This function is called by
%  o mom.run
% -------------------------------------------------------------------------
% This function calls
% o mom.optimal_weighting_matrix
% o mom.display_estimation_results_table
% o dynare_minimize_objective
% o mom.objective_function
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

if size(options_mom_.mom.weighting_matrix,1)>1 && ~(any(strcmpi('diagonal',options_mom_.mom.weighting_matrix)) || any(strcmpi('optimal',options_mom_.mom.weighting_matrix)))
    fprintf('\nYou did not specify the use of an optimal or diagonal weighting matrix. There is no point in running an iterated method of moments.\n')
end

for stage_iter = 1:size(options_mom_.mom.weighting_matrix,1)
    fprintf('Estimation stage %u\n',stage_iter);
    Woptflag = false;
    switch lower(options_mom_.mom.weighting_matrix{stage_iter})
        case 'identity_matrix'
            fprintf('  - identity weighting matrix\n');
            weighting_matrix = eye(options_mom_.mom.mom_nbr);
        case 'diagonal'
            fprintf('  - diagonal of optimal weighting matrix (Bartlett kernel with %d lags)\n', options_mom_.mom.bartlett_kernel_lag);
            if stage_iter == 1
                fprintf('    and using data-moments as initial estimate of model-moments\n');
                weighting_matrix = diag(diag(  mom.optimal_weighting_matrix(oo_.mom.m_data, oo_.mom.data_moments, options_mom_.mom.bartlett_kernel_lag)  ));
            else
                fprintf('    and using previous stage estimate of model-moments\n');
                weighting_matrix = diag(diag(  mom.optimal_weighting_matrix(oo_.mom.m_data, oo_.mom.model_moments, options_mom_.mom.bartlett_kernel_lag)  ));
            end
        case 'optimal'
            fprintf('  - optimal weighting matrix (Bartlett kernel with %d lags)\n', options_mom_.mom.bartlett_kernel_lag);
            if stage_iter == 1
                fprintf('    and using data-moments as initial estimate of model-moments\n');
                weighting_matrix = mom.optimal_weighting_matrix(oo_.mom.m_data, oo_.mom.data_moments, options_mom_.mom.bartlett_kernel_lag);
            else
                fprintf('    and using previous stage estimate of model-moments\n');
                weighting_matrix = mom.optimal_weighting_matrix(oo_.mom.m_data, oo_.mom.model_moments, options_mom_.mom.bartlett_kernel_lag);
                Woptflag = true;
            end
        otherwise % user specified matrix in file
            fprintf('  - user-specified weighting matrix\n');
            try
                load(options_mom_.mom.weighting_matrix{stage_iter},'weighting_matrix')
            catch
                error(['method_of_moments: No matrix named ''weighting_matrix'' could be found in ',options_mom_.mom.weighting_matrix{stage_iter},'.mat !'])
            end
            [nrow, ncol] = size(weighting_matrix);
            if ~isequal(nrow,ncol) || ~isequal(nrow,length(oo_.mom.data_moments)) %check if square and right size
                error(['method_of_moments: ''weighting_matrix'' must be square and have ',num2str(length(oo_.mom.data_moments)),' rows and columns!'])
            end
    end
    try % check for positive definiteness of weighting_matrix
        oo_.mom.Sw = chol(weighting_matrix);
    catch
        error('method_of_moments: Specified ''weighting_matrix'' is not positive definite. Check whether your model implies stochastic singularity!')
    end

    for optim_iter = 1:length(options_mom_.optimizer_vec)
        options_mom_.current_optimizer = options_mom_.optimizer_vec{optim_iter};
        if options_mom_.optimizer_vec{optim_iter} == 0
            xparam1 = xparam0; % no minimization, evaluate objective at current values
            fval = feval(objective_function, xparam1, BoundsInfo, oo_, estim_params_, M_, options_mom_);
        else
            if options_mom_.optimizer_vec{optim_iter} == 13
                options_mom_.mom.vector_output = true;
            else
                options_mom_.mom.vector_output = false;
            end
            if strcmp(options_mom_.mom.mom_method,'GMM') && options_mom_.mom.analytic_jacobian && ismember(options_mom_.optimizer_vec{optim_iter},options_mom_.mom.analytic_jacobian_optimizers) %do this only for gradient-based optimizers
                options_mom_.mom.compute_derivs = true;
            else
                options_mom_.mom.compute_derivs = false;
            end
            [xparam1, fval] = dynare_minimize_objective(objective_function, xparam0, options_mom_.optimizer_vec{optim_iter}, options_mom_, [BoundsInfo.lb BoundsInfo.ub], bayestopt_.name, bayestopt_, [],...
                                                                  BoundsInfo, oo_, estim_params_, M_, options_mom_);
            if options_mom_.mom.vector_output
                fval = fval'*fval;
            end
        end
        fprintf('\nStage %d Iteration %d: Value of minimized moment distance objective function: %12.10f.\n',stage_iter,optim_iter,fval)
        if options_mom_.mom.verbose
            oo_.mom = display_estimation_results_table(xparam1,NaN(size(xparam1)),M_,options_mom_,estim_params_,bayestopt_,oo_.mom,prior_dist_names,sprintf('%s (STAGE %d ITERATION %d) VERBOSE',options_mom_.mom.mom_method,stage_iter,optim_iter),sprintf('verbose_%s_stage_%d_iter_%d',lower(options_mom_.mom.mom_method),stage_iter,optim_iter));
        end
        xparam0 = xparam1;
    end
    options_mom_.vector_output = false;
    [~, ~, ~,~,~, oo_] = feval(objective_function, xparam1, BoundsInfo, oo_, estim_params_, M_, options_mom_); % get oo_.mom.model_moments for iterated GMM/SMM to compute optimal weighting matrix
end
