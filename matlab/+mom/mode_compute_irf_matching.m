function [xparam1, hessian_xparam1, fval, mom_verbose] = mode_compute_irf_matching(xparam0, hessian_xparam0, objective_function, doBayesianEstimation, weighting_info, data_moments, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% [xparam1, hessian_xparam1, fval, mom_verbose] = mode_compute_irf_matching(xparam0, hessian_xparam0, objective_function, doBayesianEstimation, weighting_info, data_moments, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% -------------------------------------------------------------------------
% Computes the minimum of the objective function (distance between data IRFs
% and model IRFs) for a sequence of optimizers.
% Note that we call a "mode" the minimum of the objective function, i.e.
% the parameter vector that minimizes the distance between the IRFs
% computed from the model and the IRFs computed from the data.
% -------------------------------------------------------------------------
% INPUTS
% xparam0:               [vector]       initialized parameters
% hessian_xparam0:       [matrix]       initialized hessian at xparam0
% objective_function:    [func handle]  name of the objective function
% doBayesianEstimation:  [logical]      true if Bayesian estimation
% weighting_info:        [structure]    information on weighting matrix
% data_moments:          [vector]       data moments
% options_mom_:          [structure]    options
% M_:                    [structure]    model information
% estim_params_:         [structure]    information on estimated parameters
% bayestopt_:            [structure]    information on priors
% BoundsInfo:            [structure]    bounds for optimization
% dr:                    [structure]    information reduced-form model
% endo_steady_state:     [vector]       steady state of endogenous variables (initval)
% exo_steady_state:      [vector]       steady state of exogenous variables (initval)
% exo_det_steady_state:  [vector]       steady state of deterministic exogenous variables (initval)
% -------------------------------------------------------------------------
% OUTPUT
% xparam1:               [vector]       mode of objective function
% hessian_xparam1:       [matrix]       hessian at xparam1
% fval:                  [double]       function value at mode
% mom_verbose:           [structure]    information on intermediate estimation results
% Also saves the computed mode and hessian to a file.
% -------------------------------------------------------------------------
% This function is called by
%  o mom.run
% -------------------------------------------------------------------------
% This function calls
%  o display_estimation_results_table
%  o dynare_minimize_objective
%  o hessian
%  o mom.objective_function
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
for optim_iter = 1:length(options_mom_.optimizer_vec)
    options_mom_.current_optimizer = options_mom_.optimizer_vec{optim_iter};
    if options_mom_.optimizer_vec{optim_iter}==0
        % no minimization, evaluate objective at current values
        xparam1 = xparam0;
        hessian_xparam1 = hessian_xparam0;
        fval = feval(objective_function, xparam1, data_moments, weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
    else
        [xparam1, fval, exitflag, hessian_xparam1, options_mom_, Scale, new_rat_hess_info] = dynare_minimize_objective(objective_function, xparam0, options_mom_.optimizer_vec{optim_iter}, options_mom_, [BoundsInfo.lb BoundsInfo.ub], bayestopt_.name, bayestopt_, hessian_xparam0,...
                                                                                                                       data_moments, weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state);        
    end
    fprintf('\nMode Compute Iteration %d: Value of minimized moment distance objective function: %12.10f.\n',optim_iter,fval);
    if options_mom_.mom.verbose
        fprintf('\n''verbose'' option: ');
        if options_mom_.cova_compute
            if options_mom_.optimizer_vec{optim_iter}==0
                hessian_xparam1_iter = hessian_xparam1;
            else
                fprintf('computing Hessian');
                hessian_xparam1_iter = hessian(objective_function, xparam1, options_mom_.gstep,...
                                               data_moments, weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
                hessian_xparam1_iter = reshape(hessian_xparam1_iter, length(xparam1), length(xparam1));
            end
            hsd_iter = sqrt(diag(hessian_xparam1_iter));
            invhessian_xparam1_iter = inv(hessian_xparam1_iter./(hsd_iter*hsd_iter'))./(hsd_iter*hsd_iter');
            std_via_invhessian_xparam1_iter = sqrt(diag(invhessian_xparam1_iter));
        else
            std_via_invhessian_xparam1_iter = NaN(size(xparam1));
        end
        fprintf(' and displaying intermediate results.');
        if doBayesianEstimation
            tbl_title_iter = sprintf('BAYESIAN %s (OPTIM ITERATION %d) VERBOSE',strrep(options_mom_.mom.mom_method,'_',' '),optim_iter);
            field_name_iter = sprintf('posterior_iter_%d',optim_iter);
        else
            tbl_title_iter = sprintf('FREQUENTIST %s (OPTIM ITERATION %d) VERBOSE',strrep(options_mom_.mom.mom_method,'_',' '),optim_iter);
            field_name_iter = sprintf('iter_%d',optim_iter);
        end
        mom_verbose.(field_name_iter) = display_estimation_results_table(xparam1,std_via_invhessian_xparam1_iter,M_,options_mom_,estim_params_,bayestopt_,[],prior_dist_names,tbl_title_iter,field_name_iter);
    end
    xparam0 = xparam1;
    hessian_xparam0 = hessian_xparam1;
end

if options_mom_.cova_compute
    if options_mom_.mom.verbose
        hessian_xparam1 = hessian_xparam1_iter;
    else
        fprintf('\nComputing Hessian at the mode.\n');
        hessian_xparam1 = hessian(objective_function, xparam1, options_mom_.gstep,...
                                  data_moments, weighting_info, options_mom_, M_, estim_params_, bayestopt_, BoundsInfo, dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
        hessian_xparam1 = reshape(hessian_xparam1, length(xparam1), length(xparam1));
    end
end
parameter_names = bayestopt_.name;
if options_mom_.cova_compute || options_mom_.mode_compute==5 || options_mom_.mode_compute==6
    hh = hessian_xparam1;
    save([M_.dname filesep 'method_of_moments' filesep M_.fname '_mode.mat'],'xparam1','hh','parameter_names','fval');
else
    save([M_.dname filesep 'method_of_moments' filesep M_.fname '_mode.mat'],'xparam1','parameter_names','fval');
end