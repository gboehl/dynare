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
%   * get_first_order_solution_params_deriv.m (previously getH.m)
%   * get_identification_jacobians.m (previously getJJ.m)
% -------------------------------------------------------------------------
% This function calls
%   * [M.fname,'.dynamic']
%   * dynare_resolve
%   * dyn_vech
%   * lyapunov_symm
%   * vec
% =========================================================================
% Copyright (C) 2011-2019 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
% =========================================================================

if nargin < 11 || isempty(useautocorr)
    useautocorr = 0;
end
if nargin < 12 || isempty(nlags)
    nlags = 3;
end
if nargin < 13 || isempty(grid_nbr)
    grid_nbr = 0;
end

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
[~,info,M,options,oo] = resol(0,M,options,oo);
options = rmfield(options,'options_ident');
oo.dr = pruned_state_space_system(M, options, oo.dr);
A = oo.dr.pruned.A;
B = oo.dr.pruned.B;
C = oo.dr.pruned.C(indvar,:);
D = oo.dr.pruned.D(indvar,:);
Om_z = oo.dr.pruned.Om_z;
Om_y = oo.dr.pruned.Om_y(indvar,indvar);
Varinov = oo.dr.pruned.Varinov;
obs_nbr = size(C,1);
%% out = [vech(cov(Y_t,Y_t)); vec(cov(Y_t,Y_{t-1}); ...; vec(cov(Y_t,Y_{t-nlags})] of indvar variables, in DR order. This is Iskrev (2010)'s J matrix.
if outputflag == 1
    % Denote Ezz0 = E_t(z_t * z_t'), then the following Lyapunov equation defines the autocovariagrom: Ezz0 -A*Ezz*A' = B*Sig_e*B'
    [Ezz0,u] = lyapunov_symm(A, Om_z, options.lyapunov_fixed_point_tol, options.qz_criterium, options.lyapunov_complex_threshold, 1, options.debug);
    stationary_vars = (1:size(C,1))';
    if ~isempty(u)
        x = abs(C*u);
        stationary_vars = find(all(x < options.Schur_vec_tol,2));
    end
    Eyy0 = NaN*ones(obs_nbr,obs_nbr);
    Eyy0(stationary_vars,stationary_vars) = C(stationary_vars,:)*Ezz0*C(stationary_vars,:)' + Om_y(stationary_vars,stationary_vars);
    indzeros = find(abs(Eyy0) < 1e-12); %find values that are numerical zero
    Eyy0(indzeros) = 0;
    if useautocorr
        sy = sqrt(diag(Ezz0)); %theoretical standard deviation
        sy = sy(stationary_vars);
        sy = sy*sy';          %cross products of standard deviations
        sy0 = sy-diag(diag(sy))+eye(length(sy));
        Eyy0corr = NaN*ones(size(C,1),size(C,1));
        Eyy0corr(stationary_vars,stationary_vars) = Eyy0./sy0;
        out = dyn_vech(Eyy0corr); %focus only on unique terms
    else
        out = dyn_vech(Eyy0); %focus only on unique terms
    end
    % compute autocovariances/autocorrelations of lagged observed variables
    tmpEyyi = A*Ezz0*C(stationary_vars,:)' + B*Varinov*D(stationary_vars,:)';
    Ai = eye(size(A,1)); %this is A^0
    for i = 1:nlags
        Eyyi = NaN*ones(obs_nbr,obs_nbr);
        Eyyi(stationary_vars,stationary_vars) = C(stationary_vars,:)*Ai*tmpEyyi;
        if useautocorr
            Eyyi = Eyyi./sy;
        end
        out = [out;vec(Eyyi)];
        Ai = Ai*A; %note that this is A^(i-1)
    end
end

%% out = vec(g_omega). This is needed for Qu and Tkachenko (2012)'s G matrix.
if outputflag == 2
% This computes the spectral density g_omega where the interval [-pi;\pi] is discretized by grid_nbr points
    freqs = (0 : pi/(grid_nbr/2):pi);% we focus only on positive values including the 0 frequency
    tpos  = exp( sqrt(-1)*freqs); %Fourier frequencies
    IA = eye(size(A,1));
    var_nbr = size(C,1);
    out = zeros(var_nbr^2*length(freqs),1);
    kk = 0;
    for ig = 1:length(freqs)
        Transferfct = D + C*((tpos(ig)*IA-A)\B);
        g_omega = (1/(2*pi))*(Transferfct*Varinov*Transferfct'); % note that ' is the conjugate transpose
        kk = kk+1;
        out(1 + (kk-1)*var_nbr^2 : kk*var_nbr^2) = g_omega(:);
    end    
end