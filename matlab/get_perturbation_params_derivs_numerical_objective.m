function [out,info] = get_perturbation_params_derivs_numerical_objective(params, outputflag, estim_params, M, oo, options)
%function [out,info] = get_perturbation_params_derivs_numerical_objective(params, outputflag, estim_params, M, oo, options)
% -------------------------------------------------------------------------
% Objective function to compute numerically the Jacobians used for get_perturbation_params_derivs
% =========================================================================
% INPUTS
%   params:         [vector] parameter values at which to evaluate objective function
%                   stderr parameters come first, corr parameters second, model parameters third
%   outputflag:     [string] flag which objective to compute (see below)
%   estim_params:   [structure] storing the estimation information
%   M:              [structure] storing the model information
%   oo:             [structure] storing the solution results
%   options:        [structure] storing the options
% -------------------------------------------------------------------------
%
% OUTPUT 
%   out (dependent on outputflag and order of approximation):
%     - 'perturbation_solution':  out = out1 = [vec(Sigma_e);vec(ghx);vec(ghu)]; (order==1)
%                                 out = out2 = [out1;vec(ghxx);vec(ghxu);vec(ghuu);vec(ghs2)]; (order==2)
%                                 out = out3 = [out1;out2;vec(ghxxx);vec(ghxxu);vec(ghxuu);vec(ghuuu);vec(ghxss);vec(ghuss)]; (order==3)
%     - 'dynamic_model':          out = [Yss; vec(g1); vec(g2); vec(g3)]
%     - 'Kalman_Transition':      out = [Yss; vec(KalmanA); dyn_vech(KalmanB*Sigma_e*KalmanB')];
%     all in DR-order
%   info            [integer] output from resol
% -------------------------------------------------------------------------
% This function is called by
%   * get_perturbation_params_derivs.m (previously getH.m)
% -------------------------------------------------------------------------
% This function calls
%   * [M.fname,'.dynamic']
%   * resol
%   * dyn_vech
% =========================================================================
% Copyright (C) 2019-2020 Dynare Team
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

%% Update stderr, corr and model parameters and compute perturbation approximation and steady state with updated parameters
M = set_all_parameters(params,estim_params,M);
[~,info,M,oo] = compute_decision_rules(M,options,oo);
Sigma_e = M.Sigma_e;

if info(1) > 0
    % there are errors in the solution algorithm
    out = [];
    return
else
    ys = oo.dr.ys; %steady state of model variables in declaration order
    ghx = oo.dr.ghx; ghu = oo.dr.ghu;
    if options.order > 1
        ghxx = oo.dr.ghxx; ghxu = oo.dr.ghxu; ghuu = oo.dr.ghuu; ghs2 = oo.dr.ghs2;
    end
    if options.order > 2
        ghxxx = oo.dr.ghxxx; ghxxu = oo.dr.ghxxu; ghxuu = oo.dr.ghxuu; ghxss = oo.dr.ghxss; ghuuu = oo.dr.ghuuu; ghuss = oo.dr.ghuss;
    end
end
Yss = ys(oo.dr.order_var); %steady state of model variables in DR order

%% out = [vec(Sigma_e);vec(ghx);vec(ghu);vec(ghxx);vec(ghxu);vec(ghuu);vec(ghs2);vec(ghxxx);vec(ghxxu);vec(ghxuu);vec(ghuuu);vec(ghxss);vec(ghuss)]
if strcmp(outputflag,'perturbation_solution')
    out = [Sigma_e(:); ghx(:); ghu(:)];
    if options.order > 1
        out = [out; ghxx(:); ghxu(:); ghuu(:); ghs2(:);];
    end
    if options.order > 2
        out = [out; ghxxx(:); ghxxu(:); ghxuu(:); ghuuu(:); ghxss(:); ghuss(:)];
    end
end

%% out = [Yss; vec(g1); vec(g2); vec(g3)]; of all endogenous variables, in DR order
if strcmp(outputflag,'dynamic_model')
    [I,~] = find(M.lead_lag_incidence'); %I is used to evaluate dynamic model files
    if options.order == 1
        [~, g1] = feval([M.fname,'.dynamic'], ys(I), oo.exo_steady_state', M.params, ys, 1);
        out = [Yss; g1(:)];
    elseif options.order == 2
        [~, g1, g2] = feval([M.fname,'.dynamic'], ys(I), oo.exo_steady_state', M.params, ys, 1);
        out = [Yss; g1(:); g2(:)];
    elseif options.order == 3
        [~, g1, g2, g3] = feval([M.fname,'.dynamic'], ys(I), oo.exo_steady_state', M.params, ys, 1);
        g3 = unfold_g3(g3, length(ys(I))+M.exo_nbr);
        out = [Yss; g1(:); g2(:); g3(:)];
    end
end

%% out = [Yss; vec(KalmanA); dyn_vech(KalmanB*Sigma_e*KalmanB')]; in DR order, where A and B are Kalman transition matrices
if strcmp(outputflag,'Kalman_Transition')
    if options.order == 1
        KalmanA = zeros(M.endo_nbr,M.endo_nbr);
        KalmanA(:,M.nstatic+(1:M.nspred)) = ghx;
        Om = ghu*Sigma_e*transpose(ghu);
        out = [Yss; KalmanA(:); dyn_vech(Om)];
    else
        error('''get_perturbation_params_derivs_numerical_objective.m'': Kalman_Transition works only at order=1');
    end
end
