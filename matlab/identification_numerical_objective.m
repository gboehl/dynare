function out = identification_numerical_objective(params, outputflag, estim_params, M, oo, options, indpmodel, indpstderr, indpcorr, indvar, useautocorr, nlags, grid_nbr)
%function out = identification_numerical_objective(params, outputflag, estim_params, M, oo, options, indpmodel, indpstderr, indpcorr, indvar, useautocorr, nlags, grid_nbr)
% -------------------------------------------------------------------------
% Objective function to compute numerically the Jacobians used for identification analysis
% Previously this function was called thet2tau.m
% =========================================================================
% INPUTS
%   params:         [vector] parameter values at which to evaluate objective function
%                   stderr parameters come first, corr parameters second, model parameters third
%   outputflag:     [integer] flag which objective to compute (see below)
%   estim_params:   [structure] storing the estimation information
%   M:              [structure] storing the model information
%   oo:             [structure] storing the reduced form solution results
%   options:        [structure] storing the options
%   indpmodel:      [vector] Index of model parameters
%   indpstderr:     [vector] Index of stderr parameters
%   indpcorr:       [matrix] Index of corr parameters
%   indvar:         [vector] Index of selected or observed variables
% -------------------------------------------------------------------------
% OUTPUTS
%   out:    dependent on outputflag
%   *  0:   out = [Yss; vec(A); vec(B); dyn_vech(Sig_e)]; of indvar variables only, in DR order. This is needed to compute dTAU and Komunjer and Ng's D.
%           Note that Jacobian of Om is computed in get_identification_Jacobians.m (previously getJJ.m) or get_first_order_solution_params_deriv.m (previously getH.m) from Jacobian of B and Sigma_e, because this is more efficient due to some testing with analytical derivatives from An and Schorfheide model
%   *  1:   out = [vech(cov(Y_t,Y_t)); vec(cov(Y_t,Y_{t-1}); ...; vec(cov(Y_t,Y_{t-nlags})] of indvar variables, in DR order. This is needed to compute Iskrev's J.
%   *  2:   out = vec(spectral density) with dimension [var_nbr^2*grid_nbr,1] Spectral density of indvar variables evaluated at (grid_nbr/2+1) discretized points in the interval [0;pi]. This is needed for Qu and Tkachenko's G.
%   * -1:   out = g1(:); of all variables, in DR order. This is needed to compute dLRE.
%   * -2:   out = [Yss; vec(A); dyn_vech(B*Sigma_e*B')]; of indvar variables only, in DR order. This is used to compute numerically second derivatives d2A, d2Om d2Yss in get_first_order_solution_params_deriv.m (previously getH.m) for kronflag=1
% where Yss is steady in DR order, A and B solution matrices of Kalman
% transition equation, Sig_e the covariance of exogenous shocks, g1 the
% Jacobian of the dynamic model equations, and Y_t selected variables
% -------------------------------------------------------------------------
% This function is called by
%   * get_identification_jacobians.m (previously getJJ.m)
% -------------------------------------------------------------------------
% This function calls
%   * [M.fname,'.dynamic']
%   * dyn_vech
%   * resol
%   * vec
% =========================================================================
% Copyright (C) 2011-2020 Dynare Team
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

%% Update stderr, corr and model parameters
%note that if no estimated_params_block is given, then all stderr and model parameters are selected but no corr parameters
if length(params) > length(indpmodel)
    if isempty(indpstderr)==0 && isempty(estim_params.var_exo) %if there are stderr parameters but no estimated_params_block
        %provide temporary necessary information for stderr parameters
        estim_params.nvx = length(indpstderr);
        estim_params.var_exo = indpstderr';
    end
    if isempty(indpmodel)==0 && isempty(estim_params.param_vals) %if there are model parameters but no estimated_params_block
        %provide temporary necessary information for model parameters
        estim_params.np = length(indpmodel);
        estim_params.param_vals = indpmodel';
    end
    M = set_all_parameters(params,estim_params,M); %this function can only be used if there is some information in estim_params
else
    %if there are only model parameters, we don't need to use set_all_parameters
    M.params(indpmodel) = params;
end

%% compute Kalman transition matrices and steady state with updated parameters
[~,info,M,oo] = compute_decision_rules(M,options,oo);
options = rmfield(options,'options_ident');
pruned = pruned_state_space_system(M, options, oo.dr, indvar, nlags, useautocorr, 0);

%% out = [vech(cov(Y_t,Y_t)); vec(cov(Y_t,Y_{t-1}); ...; vec(cov(Y_t,Y_{t-nlags})] of indvar variables, in DR order. This is Iskrev (2010)'s J matrix.
if outputflag == 1    
    out = dyn_vech(pruned.Var_y);
    for i = 1:nlags
        if useautocorr
            out = [out;vec(pruned.Corr_yi(:,:,i))];
        else
            out = [out;vec(pruned.Var_yi(:,:,i))];
        end
    end
end

%% out = vec(g_omega). This is needed for Qu and Tkachenko (2012)'s G matrix.
if outputflag == 2
    % This computes the spectral density g_omega where the interval [-pi;\pi] is discretized by grid_nbr points
    freqs = (0 : pi/(grid_nbr/2):pi);% we focus only on positive values including the 0 frequency
    tpos  = exp( sqrt(-1)*freqs); %Fourier frequencies
    IA = eye(size(pruned.A,1));
    var_nbr = size(pruned.C,1);
    out = zeros(var_nbr^2*length(freqs),1);
    kk = 0;
    for ig = 1:length(freqs)
        Transferfct = pruned.D + pruned.C*((tpos(ig)*IA-pruned.A)\pruned.B);
        g_omega = (1/(2*pi))*(Transferfct*pruned.Varinov*Transferfct'); % note that ' is the conjugate transpose
        kk = kk+1;
        out(1 + (kk-1)*var_nbr^2 : kk*var_nbr^2) = g_omega(:);
    end
end