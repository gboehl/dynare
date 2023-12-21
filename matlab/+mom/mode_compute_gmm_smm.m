function [xparam1, weighting_info, mom_verbose] = mode_compute_gmm_smm(xparam0, objective_function, m_data, data_moments, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% [xparam1, weighting_info, mom_verbose] = mode_compute_gmm_smm(xparam0, objective_function, m_data, data_moments, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% -------------------------------------------------------------------------
% Iterated method of moments for GMM and SMM, computes the minimum of the
% objective function (distance between data moments and model moments)
% for a sequence of optimizers and GMM/SMM iterations with different
% weighting matrices.
% -------------------------------------------------------------------------
% INPUTS
% xparam0:               [vector]       vector of initialized parameters
% objective_function:    [func handle]  name of the objective function
% m_data:                [matrix]       selected data moments at each point in time
% data_moments:          [vector]       vector of data moments
% options_mom_:          [structure]    options
% M_:                    [structure]    model information
% estim_params_:         [structure]    information on estimated parameters
% bayestopt_:            [structure]    information on priors
% BoundsInfo:            [structure]    bounds for optimization
% dr:                    [structure]    reduced form model
% endo_steady_state:     [vector]       steady state for endogenous variables (initval)
% exo_steady_state:      [vector]       steady state for exogenous variables (initval)
% exo_det_steady_state:  [vector]       steady state for exogenous deterministic variables (initval)
% -------------------------------------------------------------------------
% OUTPUT
% xparam1:               [vector]       mode of objective function
% weighting_info:        [structure]    information on weighting matrix
% mom_verbose:           [structure]    information on intermediate estimation results
% -------------------------------------------------------------------------
% This function is called by
%  o mom.run
% -------------------------------------------------------------------------
% This function calls
% o mom.optimal_weighting_matrix
% o mom.display_estimation_results_table
% o dynare_minimize_objective
% o mom.objective_function
% o prior_dist_names
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


mom_verbose = [];
if size(options_mom_.mom.weighting_matrix,1)>1 && ~(any(strcmpi('diagonal',options_mom_.mom.weighting_matrix)) || any(strcmpi('optimal',options_mom_.mom.weighting_matrix)))
    fprintf('\nYou did not specify the use of an optimal or diagonal weighting matrix. There is no point in running an iterated method of moments.\n');
end

for stage_iter = 1:size(options_mom_.mom.weighting_matrix,1)
    fprintf('Estimation stage %u\n',stage_iter);
    weighting_info.Woptflag = false;
    switch lower(options_mom_.mom.weighting_matrix{stage_iter})
        case 'identity_matrix'
            fprintf('  - identity weighting matrix\n');
            weighting_matrix = eye(options_mom_.mom.mom_nbr);
        case 'diagonal'
            fprintf('  - diagonal of optimal weighting matrix (Bartlett kernel with %d lags)\n', options_mom_.mom.bartlett_kernel_lag);
            if stage_iter == 1
                fprintf('    and using data-moments as initial estimate of model-moments\n');
                weighting_matrix = diag(diag(  mom.optimal_weighting_matrix(m_data, data_moments, options_mom_.mom.bartlett_kernel_lag)  ));
            else
                fprintf('    and using previous stage estimate of model-moments\n');
                weighting_matrix = diag(diag(  mom.optimal_weighting_matrix(m_data, model_moments, options_mom_.mom.bartlett_kernel_lag)  ));
            end
        case 'optimal'
            fprintf('  - optimal weighting matrix (Bartlett kernel with %d lags)\n', options_mom_.mom.bartlett_kernel_lag);
            if stage_iter == 1
                fprintf('    and using data-moments as initial estimate of model-moments\n');
                weighting_matrix = mom.optimal_weighting_matrix(m_data, data_moments, options_mom_.mom.bartlett_kernel_lag);
            else
                fprintf('    and using previous stage estimate of model-moments\n');
                weighting_matrix = mom.optimal_weighting_matrix(m_data, model_moments, options_mom_.mom.bartlett_kernel_lag);
                weighting_info.Woptflag = true;
            end
        otherwise % user specified matrix in file
            fprintf('  - user-specified weighting matrix\n');
            try
                load(options_mom_.mom.weighting_matrix{stage_iter},'weighting_matrix')
            catch
                error(['method_of_moments: No matrix named ''weighting_matrix'' could be found in ',options_mom_.mom.weighting_matrix{stage_iter},'.mat !']);
            end
            [nrow, ncol] = size(weighting_matrix);
            if ~isequal(nrow,ncol) || ~isequal(nrow,length(data_moments)) %check if square and right size
                error(['method_of_moments: ''weighting_matrix'' must be square and have ',num2str(length(data_moments)),' rows and columns!']);
            end
    end
    try % check for positive definiteness of weighting_matrix
        weighting_info.Sw = chol(weighting_matrix);
    catch
        error('method_of_moments: Specified ''weighting_matrix'' is not positive definite. Check whether your model implies stochastic singularity!');
    end

    for optim_iter = 1:length(options_mom_.optimizer_vec)
        options_mom_.current_optimizer = options_mom_.optimizer_vec{optim_iter};
        if options_mom_.optimizer_vec{optim_iter} == 0
            xparam1 = xparam0; % no minimization, evaluate objective at current values            
            fval = feval(objective_function, xparam1, data_moments, weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
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
                                                                  data_moments, weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
            if options_mom_.mom.vector_output
                fval = fval'*fval;
            end
        end
        fprintf('\nStage %d Iteration %d: Value of minimized moment distance objective function: %12.10f.\n',stage_iter,optim_iter,fval);
        if options_mom_.mom.verbose
            fprintf('\n''verbose'' option: ');
            std_via_invhessian_xparam1_iter = NaN(size(xparam1));            
            tbl_title_iter = sprintf('FREQUENTIST %s (STAGE %d ITERATION %d) VERBOSE',options_mom_.mom.mom_method,stage_iter,optim_iter);
            field_name_iter = sprintf('%s_stage_%d_iter_%d',lower(options_mom_.mom.mom_method),stage_iter,optim_iter);
            mom_verbose.(field_name_iter) = display_estimation_results_table(xparam1,std_via_invhessian_xparam1_iter,M_,options_mom_,estim_params_,bayestopt_,[],prior_dist_names,tbl_title_iter,field_name_iter);
        end
        xparam0 = xparam1;
    end
    options_mom_.vector_output = false;
    [~, ~, ~, ~, ~, ~, model_moments] = feval(objective_function, xparam1, data_moments, weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state); % get model_moments for iterated GMM/SMM to compute optimal weighting matrix
end
