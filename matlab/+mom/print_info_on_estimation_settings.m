function print_info_on_estimation_settings(options_mom_, number_of_estimated_parameters, do_bayesian_estimation)
% print_info_on_estimation_settings(options_mom_, number_of_estimated_parameters, do_bayesian_estimation)
% -------------------------------------------------------------------------
% Print information on the method of moments estimation settings to the console
% -------------------------------------------------------------------------
% INPUTS
% options_mom_                    [struct]   options for the method of moments estimation
% number_of_estimated_parameters  [integer]  number of estimated parameters
% do_bayesian_estimation          [boolean]  true if the estimation is Bayesian
% -------------------------------------------------------------------------
% OUTPUT
% No output, just displays the chosen settings
% -------------------------------------------------------------------------
% This function is called by
%  o mom.run
% -------------------------------------------------------------------------
% This function calls
%  o skipline
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


fprintf('\n---------------------------------------------------\n')
if strcmp(options_mom_.mom.mom_method,'SMM')
    fprintf('Simulated method of moments with');
elseif strcmp(options_mom_.mom.mom_method,'GMM')
    fprintf('General method of moments with');
end
if strcmp(options_mom_.mom.mom_method,'SMM') || strcmp(options_mom_.mom.mom_method,'GMM')
    if options_mom_.prefilter
        fprintf('\n  - centered moments (prefilter=1)');
    else
        fprintf('\n  - uncentered moments (prefilter=0)');
    end
    if options_mom_.mom.penalized_estimator
        fprintf('\n  - penalized estimation using deviation from prior mean and weighted with prior precision');
    end
end
if strcmp(options_mom_.mom.mom_method,'IRF_MATCHING')
    if do_bayesian_estimation
        fprintf('Bayesian Impulse Response Function Matching with');
    else
        fprintf('Frequentist Impulse Response Function Matching with');
    end
    if ~isempty(options_mom_.mom.irf_matching_file.name)
        fprintf('\n  - irf_matching_file: %s',[options_mom_.mom.irf_matching_file.path filesep options_mom_.mom.irf_matching_file.name '.m']);
    end    
end
for i = 1:length(options_mom_.optimizer_vec)
    if i == 1
        str = '- optimizer (mode_compute';
    else
        str = '            (additional_optimizer_steps';
    end
    switch options_mom_.optimizer_vec{i}
        case 0
            fprintf('\n  %s=0): no minimization',str);
        case 1
            fprintf('\n  %s=1): fmincon',str);
        case 2
            fprintf('\n  %s=2): continuous simulated annealing',str);
        case 3
            fprintf('\n  %s=3): fminunc',str);
        case 4
            fprintf('\n  %s=4): csminwel',str);
        case 5
            fprintf('\n  %s=5): newrat',str);
        case 6
            fprintf('\n  %s=6): gmhmaxlik',str);
        case 7
            fprintf('\n  %s=7): fminsearch',str);
        case 8
            fprintf('\n  %s=8): Dynare Nelder-Mead simplex',str);
        case 9
            fprintf('\n  %s=9): CMA-ES',str);
        case 10
            fprintf('\n  %s=10): simpsa',str);
        case 11
            skipline;
            error('method_of_moments: online_auxiliary_filter (mode_compute=11) is only supported with likelihood-based estimation techniques!');
        case 12
            fprintf('\n  %s=12): particleswarm',str);
        case 101
            fprintf('\n  %s=101): SolveOpt',str);
        case 102
            fprintf('\n  %s=102): simulannealbnd',str);
        case 13
            fprintf('\n  %s=13): lsqnonlin',str);
        otherwise
            if ischar(options_mom_.optimizer_vec{i})
                fprintf('\n  %s=%s): user-defined',str,options_mom_.optimizer_vec{i});
            else
                skipline;
                error('method_of_moments: Unknown optimizer!');
            end
    end
    if options_mom_.silent_optimizer
        fprintf(' (silent)');
    end
    if strcmp(options_mom_.mom.mom_method,'GMM') && options_mom_.mom.analytic_jacobian && ismember(options_mom_.optimizer_vec{i},options_mom_.mom.analytic_jacobian_optimizers)
        fprintf(' (using analytical Jacobian)');
    end
end
if options_mom_.order > 0
    fprintf('\n  - stochastic simulations with perturbation order: %d', options_mom_.order)
end
if options_mom_.order > 1 && options_mom_.pruning
    fprintf(' (with pruning)')
end
if strcmp(options_mom_.mom.mom_method,'GMM') || strcmp(options_mom_.mom.mom_method,'SMM')
    if strcmp(options_mom_.mom.mom_method,'GMM') && options_mom_.mom.analytic_standard_errors
        fprintf('\n  - standard errors: analytic derivatives');
    else
        fprintf('\n  - standard errors: numerical derivatives');
    end
    fprintf('\n  - number of matched moments: %d', options_mom_.mom.mom_nbr);
elseif strcmp(options_mom_.mom.mom_method,'IRF_MATCHING')
    fprintf('\n  - number of matched IRFs: %d', options_mom_.mom.mom_nbr);
end
fprintf('\n  - number of parameters: %d', number_of_estimated_parameters);
fprintf('\n\n');