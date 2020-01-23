function [MEAN, dMEAN, REDUCEDFORM, dREDUCEDFORM, DYNAMIC, dDYNAMIC, MOMENTS, dMOMENTS, dSPECTRUM, dMINIMAL, derivatives_info] = get_identification_jacobians(estim_params, M, oo, options, options_ident, indpmodel, indpstderr, indpcorr, indvobs)
% function [MEAN, dMEAN, REDUCEDFORM, dREDUCEDFORM, DYNAMIC, dDYNAMIC, MOMENTS, dMOMENTS, dSPECTRUM, dMINIMAL, derivatives_info] = get_identification_jacobians(estim_params, M, oo, options, options_ident, indpmodel, indpstderr, indpcorr, indvobs)
% previously getJJ.m in Dynare 4.5
% Sets up the Jacobians needed for identification analysis
% =========================================================================
% INPUTS
%   estim_params:   [structure] storing the estimation information
%   M:              [structure] storing the model information
%   oo:             [structure] storing the reduced-form solution results
%   options:        [structure] storing the options
%   options_ident:  [structure] storing the options for identification
%   indpmodel:      [modparam_nbr by 1] index of estimated parameters in M_.params;
%                     corresponds to model parameters (no stderr and no corr)
%                     in estimated_params block; if estimated_params block is
%                     not available, then all model parameters are selected
%   indpstderr:     [stderrparam_nbr by 1] index of estimated standard errors,
%                     i.e. for all exogenous variables where "stderr" is given
%                     in the estimated_params block; if estimated_params block
%                     is not available, then all stderr parameters are selected
%   indpcorr:       [corrparam_nbr by 2] matrix of estimated correlations,
%                     i.e. for all exogenous variables where "corr" is given
%                     in the estimated_params block; if estimated_params block
%                     is not available, then no corr parameters are selected
%   indvobs:        [obs_nbr by 1] index of observed (VAROBS) variables
% -------------------------------------------------------------------------
% OUTPUTS
%
%  MEAN             [endo_nbr by 1], in DR order. Expectation of all model variables
%                   * order==1: corresponds to steady state
%                   * order==2|3: corresponds to mean computed from pruned state space system (as in Andreasen, Fernandez-Villaverde, Rubio-Ramirez, 2018)
%  dMEAN            [endo_nbr by totparam_nbr], in DR Order, Jacobian (wrt all params) of MEAN
%
%  REDUCEDFORM      [rowredform_nbr by 1] in DR order. Steady state and reduced-form model solution matrices for all model variables
%                   * order==1: [Yss' vec(ghx)' vech(ghu*Sigma_e*ghu')']',
%                     where rowredform_nbr = endo_nbr*(1+nspred+(endo_nbr+1)/2)
%                   * order==2: [Yss' vec(ghx)' vech(ghu*Sigma_e*ghu')' vec(ghxx)' vec(ghxu)' vec(ghuu)' vec(ghs2)']',
%                     where rowredform_nbr = endo_nbr*(1+nspred+(endo_nbr+1)/2+nspred^2+nspred*exo_nr+exo_nbr^2+1)
%                   * order==3: [Yss' vec(ghx)' vech(ghu*Sigma_e*ghu')' vec(ghxx)' vec(ghxu)' vec(ghuu)' vec(ghs2)' vec(ghxxx)' vec(ghxxu)' vec(ghxuu)' vec(ghuuu)' vec(ghxss)' vec(ghuss)']',
%                     where rowredform_nbr = endo_nbr*(1+nspred+(endo_nbr+1)/2+nspred^2+nspred*exo_nr+exo_nbr^2+1+nspred^3+nspred^2*exo_nbr+nspred*exo_nbr^2+exo_nbr^3+nspred+exo_nbr)
%  dREDUCEDFORM:    [rowredform_nbr by totparam_nbr] in DR order, Jacobian (wrt all params) of REDUCEDFORM
%                   * order==1: corresponds to Iskrev (2010)'s J_2 matrix
%                   * order==2: corresponds to Mutschler (2015)'s J matrix
%
%  DYNAMIC          [rowdyn_nbr by 1] in declaration order. Steady state and dynamic model derivatives for all model variables
%                   * order==1: [ys' vec(g1)']', rowdyn_nbr=endo_nbr+length(g1)
%                   * order==2: [ys' vec(g1)' vec(g2)']', rowdyn_nbr=endo_nbr+length(g1)+length(g2)
%                   * order==3: [ys' vec(g1)' vec(g2)' vec(g3)']', rowdyn_nbr=endo_nbr+length(g1)+length(g2)+length(g3)
%  dDYNAMIC         [rowdyn_nbr by modparam_nbr] in declaration order. Jacobian (wrt model parameters) of DYNAMIC
%                   * order==1: corresponds to Ratto and Iskrev (2011)'s J_\Gamma matrix (or LRE)
%
%  MOMENTS:         [obs_nbr+obs_nbr*(obs_nbr+1)/2+nlags*obs_nbr^2 by 1] in DR order. First two theoretical moments for VAROBS variables, i.e.
%                   [E[varobs]' vech(E[varobs*varobs'])' vec(E[varobs*varobs(-1)'])' ... vec(E[varobs*varobs(-nlag)'])']
%  dMOMENTS:        [obs_nbr+obs_nbr*(obs_nbr+1)/2+nlags*obs_nbr^2 by totparam_nbr] in DR order. Jacobian (wrt all params) of MOMENTS
%                   * order==1: corresponds to Iskrev (2010)'s J matrix
%                   * order==2: corresponds to Mutschler (2015)'s \bar{M}_2 matrix, i.e. theoretical moments from the pruned state space system
%
%  dSPECTRUM:       [totparam_nbr by totparam_nbr] in DR order. Gram matrix of Jacobian (wrt all params) of spectral density for VAROBS variables, where
%                   spectral density at frequency w: f(w) = (2*pi)^(-1)*H(exp(-i*w))*E[Inov*Inov']*ctranspose(H(exp(-i*w)) with H being the Transfer function
%                   dSPECTRUM = int_{-\pi}^\pi transpose(df(w)/dp')*(df(w)/dp') dw
%                   * order==1: corresponds to Qu and Tkachenko (2012)'s G matrix, where Inov and H are computed from linear state space system
%                   * order==2: corresponds to Mutschler (2015)'s G_2 matrix, where Inov and H are computed from second-order pruned state space system
%                   * order==3: Inov and H are computed from third-order pruned state space system
%
%  dMINIMAL:        [obs_nbr+minx_nbr^2+minx_nbr*exo_nbr+obs_nbr*minx_nbr+obs_nbr*exo_nbr+exo_nbr*(exo_nbr+1)/2 by totparam_nbr+minx_nbr^2+exo_nbr^2]
%                   Jacobian (wrt all params, and similarity_transformation_matrices (T and U)) of observational equivalent minimal ABCD system,
%                   corresponds to Komunjer and Ng (2011)'s Deltabar matrix, where
%                   MINIMAL = [vec(E[varobs]' vec(minA)' vec(minB)' vec(minC)' vec(minD)' vech(Sigma_e)']'
%                   minA, minB, minC and minD is the minimal state space system computed in get_minimal_state_representation
%                   * order==1: E[varobs] is equal to steady state
%                   * order==2|3: E[varobs] is computed from the pruned state space system (second|third-order accurate), as noted in section 5 of Komunjer and Ng (2011)
%
%  derivatives_info [structure] for use in dsge_likelihood to compute Hessian analytically. Only used at order==1.
%                   Contains dA, dB, and d(B*Sigma_e*B'), where A and B are Kalman filter transition matrice.
%
% -------------------------------------------------------------------------
% This function is called by
%   * identification_analysis.m
% -------------------------------------------------------------------------
% This function calls
%   * commutation
%   * get_minimal_state_representation
%   * duplication
%   * dyn_vech
%   * get_perturbation_params_derivs (previously getH)
%   * get_all_parameters
%   * fjaco
%   * lyapunov_symm
%   * identification_numerical_objective (previously thet2tau)
%   * vec
% =========================================================================
% Copyright (C) 2010-2020 Dynare Team
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

%get fields from options_ident
no_identification_moments   = options_ident.no_identification_moments;
no_identification_minimal   = options_ident.no_identification_minimal;
no_identification_spectrum  = options_ident.no_identification_spectrum;
order       = options_ident.order;
nlags       = options_ident.ar;
useautocorr = options_ident.useautocorr;
grid_nbr    = options_ident.grid_nbr;
kronflag    = options_ident.analytic_derivation_mode;

% get fields from M
endo_nbr           = M.endo_nbr;
exo_nbr            = M.exo_nbr;
fname              = M.fname;
lead_lag_incidence = M.lead_lag_incidence;
nspred             = M.nspred;
nstatic            = M.nstatic;
params             = M.params;
Sigma_e            = M.Sigma_e;
stderr_e           = sqrt(diag(Sigma_e));

% set all selected values: stderr and corr come first, then model parameters
xparam1 = get_all_parameters(estim_params, M); %try using estimated_params block
if isempty(xparam1)
    %if there is no estimated_params block, consider all stderr and all model parameters, but no corr parameters
    xparam1 = [stderr_e', params'];
end

%get numbers/lengths of vectors
modparam_nbr    = length(indpmodel);
stderrparam_nbr = length(indpstderr);
corrparam_nbr   = size(indpcorr,1);
totparam_nbr    = stderrparam_nbr + corrparam_nbr + modparam_nbr;
obs_nbr         = length(indvobs);
d2flag          = 0; % do not compute second parameter derivatives

% Get Jacobians (wrt selected params) of steady state, dynamic model derivatives and perturbation solution matrices for all endogenous variables
oo.dr.derivs = get_perturbation_params_derivs(M, options, estim_params, oo, indpmodel, indpstderr, indpcorr, d2flag);

[I,~] = find(lead_lag_incidence'); %I is used to select nonzero columns of the Jacobian of endogenous variables in dynamic model files
yy0 = oo.dr.ys(I);           %steady state of dynamic (endogenous and auxiliary variables) in lead_lag_incidence order
Yss = oo.dr.ys(oo.dr.order_var); % steady state in DR order
if order == 1
    [~, g1 ] = feval([fname,'.dynamic'], yy0, oo.exo_steady_state', params, oo.dr.ys, 1);
    %g1 is [endo_nbr by yy0ex0_nbr first derivative (wrt all dynamic variables) of dynamic model equations, i.e. df/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
    DYNAMIC = [Yss;
               vec(g1(oo.dr.order_var,:))]; %add steady state and put rows of g1 in DR order
    dDYNAMIC = [oo.dr.derivs.dYss;
                reshape(oo.dr.derivs.dg1(oo.dr.order_var,:,:),size(oo.dr.derivs.dg1,1)*size(oo.dr.derivs.dg1,2),size(oo.dr.derivs.dg1,3)) ]; %reshape dg1 in DR order and add steady state
    REDUCEDFORM = [Yss;
                   vec(oo.dr.ghx);
                   dyn_vech(oo.dr.ghu*Sigma_e*transpose(oo.dr.ghu))]; %in DR order
    dREDUCEDFORM = zeros(endo_nbr*nspred+endo_nbr*(endo_nbr+1)/2, totparam_nbr);
    for j=1:totparam_nbr
        dREDUCEDFORM(:,j) = [vec(oo.dr.derivs.dghx(:,:,j));
                            dyn_vech(oo.dr.derivs.dOm(:,:,j))];
    end
    dREDUCEDFORM = [ [zeros(endo_nbr, stderrparam_nbr+corrparam_nbr) oo.dr.derivs.dYss]; dREDUCEDFORM ]; % add steady state

elseif order == 2
    [~, g1, g2 ] = feval([fname,'.dynamic'], yy0, oo.exo_steady_state', params, oo.dr.ys, 1);
    %g1 is [endo_nbr by yy0ex0_nbr first derivative (wrt all dynamic variables) of dynamic model equations, i.e. df/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
    %g2 is [endo_nbr by yy0ex0_nbr^2] second derivative (wrt all dynamic variables) of dynamic model equations, i.e. d(df/dyy0ex0)/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
    DYNAMIC = [Yss;
               vec(g1(oo.dr.order_var,:));
               vec(g2(oo.dr.order_var,:))]; %add steady state and put rows of g1 and g2 in DR order
    dDYNAMIC = [oo.dr.derivs.dYss;
                reshape(oo.dr.derivs.dg1(oo.dr.order_var,:,:),size(oo.dr.derivs.dg1,1)*size(oo.dr.derivs.dg1,2),size(oo.dr.derivs.dg1,3));  %reshape dg1 in DR order
                reshape(oo.dr.derivs.dg2(oo.dr.order_var,:),size(oo.dr.derivs.dg1,1)*size(oo.dr.derivs.dg1,2)^2,size(oo.dr.derivs.dg1,3))]; %reshape dg2 in DR order
    REDUCEDFORM = [Yss;
                   vec(oo.dr.ghx);
                   dyn_vech(oo.dr.ghu*Sigma_e*transpose(oo.dr.ghu));
                   vec(oo.dr.ghxx);
                   vec(oo.dr.ghxu);
                   vec(oo.dr.ghuu);
                   vec(oo.dr.ghs2)]; %in DR order
    dREDUCEDFORM = zeros(endo_nbr*nspred+endo_nbr*(endo_nbr+1)/2+endo_nbr*nspred^2+endo_nbr*nspred*exo_nbr+endo_nbr*exo_nbr^2+endo_nbr, totparam_nbr);
    for j=1:totparam_nbr
        dREDUCEDFORM(:,j) = [vec(oo.dr.derivs.dghx(:,:,j));
                            dyn_vech(oo.dr.derivs.dOm(:,:,j));
                            vec(oo.dr.derivs.dghxx(:,:,j));
                            vec(oo.dr.derivs.dghxu(:,:,j));
                            vec(oo.dr.derivs.dghuu(:,:,j));
                            vec(oo.dr.derivs.dghs2(:,j))];
    end
    dREDUCEDFORM = [ [zeros(endo_nbr, stderrparam_nbr+corrparam_nbr) oo.dr.derivs.dYss]; dREDUCEDFORM ]; % add steady state
elseif order == 3
    [~, g1, g2, g3 ] = feval([fname,'.dynamic'], yy0, oo.exo_steady_state', params, oo.dr.ys, 1);
    %g1 is [endo_nbr by yy0ex0_nbr first derivative (wrt all dynamic variables) of dynamic model equations, i.e. df/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
    %g2 is [endo_nbr by yy0ex0_nbr^2] second derivative (wrt all dynamic variables) of dynamic model equations, i.e. d(df/dyy0ex0)/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
    DYNAMIC = [Yss;
               vec(g1(oo.dr.order_var,:));
               vec(g2(oo.dr.order_var,:));
               vec(g3(oo.dr.order_var,:))]; %add steady state and put rows of g1 and g2 in DR order
    dDYNAMIC = [oo.dr.derivs.dYss;
                reshape(oo.dr.derivs.dg1(oo.dr.order_var,:,:),size(oo.dr.derivs.dg1,1)*size(oo.dr.derivs.dg1,2),size(oo.dr.derivs.dg1,3));  %reshape dg1 in DR order
                reshape(oo.dr.derivs.dg2(oo.dr.order_var,:),size(oo.dr.derivs.dg1,1)*size(oo.dr.derivs.dg1,2)^2,size(oo.dr.derivs.dg1,3));
                reshape(oo.dr.derivs.dg2(oo.dr.order_var,:),size(oo.dr.derivs.dg1,1)*size(oo.dr.derivs.dg1,2)^2,size(oo.dr.derivs.dg1,3))]; %reshape dg3 in DR order
    REDUCEDFORM = [Yss;
                   vec(oo.dr.ghx);
                   dyn_vech(oo.dr.ghu*Sigma_e*transpose(oo.dr.ghu));
                   vec(oo.dr.ghxx); vec(oo.dr.ghxu); vec(oo.dr.ghuu); vec(oo.dr.ghs2);
                   vec(oo.dr.ghxxx); vec(oo.dr.ghxxu); vec(oo.dr.ghxuu); vec(oo.dr.ghuuu); vec(oo.dr.ghxss); vec(oo.dr.ghuss)]; %in DR order
    dREDUCEDFORM = zeros(size(REDUCEDFORM,1)-endo_nbr, totparam_nbr);
    for j=1:totparam_nbr
        dREDUCEDFORM(:,j) = [vec(oo.dr.derivs.dghx(:,:,j));
                             dyn_vech(oo.dr.derivs.dOm(:,:,j));
                             vec(oo.dr.derivs.dghxx(:,:,j)); vec(oo.dr.derivs.dghxu(:,:,j)); vec(oo.dr.derivs.dghuu(:,:,j)); vec(oo.dr.derivs.dghs2(:,j))
                             vec(oo.dr.derivs.dghxxx(:,:,j)); vec(oo.dr.derivs.dghxxu(:,:,j)); vec(oo.dr.derivs.dghxuu(:,:,j)); vec(oo.dr.derivs.dghuuu(:,:,j)); vec(oo.dr.derivs.dghxss(:,:,j)); vec(oo.dr.derivs.dghuss(:,:,j))];
    end
    dREDUCEDFORM = [ [zeros(endo_nbr, stderrparam_nbr+corrparam_nbr) oo.dr.derivs.dYss]; dREDUCEDFORM ]; % add steady state
end

% Get (pruned) state space representation:
options.options_ident.indvobs = indvobs;
options.options_ident.indpmodel = indpmodel;
options.options_ident.indpstderr = indpstderr;
options.options_ident.indpcorr = indpcorr;
oo.dr = pruned_state_space_system(M, options, oo.dr);
MEAN = oo.dr.pruned.E_y;         dMEAN = oo.dr.pruned.dE_y;
A = oo.dr.pruned.A;             dA   = oo.dr.pruned.dA;
B = oo.dr.pruned.B;             dB   = oo.dr.pruned.dB;
C = oo.dr.pruned.C;             dC   = oo.dr.pruned.dC;
D = oo.dr.pruned.D;             dD   = oo.dr.pruned.dD;
c = oo.dr.pruned.c;             dc   = oo.dr.pruned.dc;
d = oo.dr.pruned.d;             dd   = oo.dr.pruned.dd;
Varinov = oo.dr.pruned.Varinov; dVarinov = oo.dr.pruned.dVarinov;
Om_z = oo.dr.pruned.Om_z; dOm_z = oo.dr.pruned.dOm_z;
Om_y = oo.dr.pruned.Om_y; dOm_y = oo.dr.pruned.dOm_y;

%storage for Jacobians used in dsge_likelihood.m for analytical Gradient and Hession of likelihood (only at order=1)
derivatives_info = struct();
if order == 1
    dT = zeros(endo_nbr,endo_nbr,totparam_nbr);
    dT(:,(nstatic+1):(nstatic+nspred),:) = oo.dr.derivs.dghx;
    derivatives_info.DT   = dT;
    derivatives_info.DOm  = oo.dr.derivs.dOm;
    derivatives_info.DYss = oo.dr.derivs.dYss;
end

%% Compute dMOMENTS
%Note that our state space system for computing second moments is the following:
% zhat = A*zhat(-1) + B*xi, where zhat = z - E(z)
% yhat = C*zhat(-1) + D*xi, where yhat = y - E(y)
if ~no_identification_moments
    MOMENTS = identification_numerical_objective(xparam1, 1, estim_params, M, oo, options, indpmodel, indpstderr, indpcorr, indvobs, useautocorr, nlags, grid_nbr); %[outputflag=1]
    MOMENTS = [MEAN; MOMENTS];
    if kronflag == -1
        %numerical derivative of autocovariogram
        dMOMENTS = fjaco(str2func('identification_numerical_objective'), xparam1, 1, estim_params, M, oo, options, indpmodel, indpstderr, indpcorr, indvobs, useautocorr, nlags, grid_nbr); %[outputflag=1]
        dMOMENTS = [dMEAN; dMOMENTS]; %add Jacobian of steady state of VAROBS variables
    else
        dMOMENTS = zeros(obs_nbr + obs_nbr*(obs_nbr+1)/2 + nlags*obs_nbr^2 , totparam_nbr);
        dMOMENTS(1:obs_nbr,:) = dMEAN; %add Jacobian of first moments of VAROBS variables
        % Denote Ezz0 = E[zhat*zhat'], then the following Lyapunov equation defines the autocovariagram: Ezz0 -A*Ezz0*A' = B*Sig_xi*B' = Om_z
        [Ezz0,u] = lyapunov_symm(A, Om_z, options.lyapunov_fixed_point_tol, options.qz_criterium, options.lyapunov_complex_threshold, 1, options.debug);
        stationary_vars = (1:length(indvobs))';
        if ~isempty(u)
            x = abs(C*u);
            stationary_vars = find(all(x < options.Schur_vec_tol,2));
        end
        Eyy0 = NaN*ones(obs_nbr,obs_nbr);
        Eyy0(stationary_vars,stationary_vars) = C(stationary_vars,:)*Ezz0*C(stationary_vars,:)' + Om_y(stationary_vars,stationary_vars);
        %here method=1 is used, whereas all other calls of lyapunov_symm use method=2. The reason is that T and U are persistent, and input matrix A is the same, so using option 2 for all the rest of iterations spares a lot of computing time while not repeating Schur every time
        indzeros = find(abs(Eyy0) < 1e-12); %find values that are numerical zero
        Eyy0(indzeros) = 0;
        if useautocorr
            sdy = sqrt(diag(Eyy0)); %theoretical standard deviation
            sdy = sdy(stationary_vars);
            sy = sdy*sdy';          %cross products of standard deviations
        end

        for jp = 1:totparam_nbr
            if jp <= (stderrparam_nbr+corrparam_nbr)
                %Note that for stderr and corr parameters, the derivatives of the system matrices are zero, i.e. dA=dB=dC=dD=0
                dEzz0 = lyapunov_symm(A,dOm_z(:,:,jp),options.lyapunov_fixed_point_tol,options.qz_criterium,options.lyapunov_complex_threshold,2,options.debug);     %here method=2 is used to spare a lot of computing time while not repeating Schur every time
                dEyy0 = C*dEzz0*C' + dOm_y(:,:,jp);
            else %model parameters
                dEzz0 = lyapunov_symm(A,dOm_z(:,:,jp)+dA(:,:,jp)*Ezz0*A'+A*Ezz0*dA(:,:,jp)',options.lyapunov_fixed_point_tol,options.qz_criterium,options.lyapunov_complex_threshold,2,options.debug);     %here method=2 is used to spare a lot of computing time while not repeating Schur every time
                dEyy0 = dC(:,:,jp)*Ezz0*C' + C*dEzz0*C' + C*Ezz0*dC(:,:,jp)' + dOm_y(:,:,jp);
            end
            %indzeros = find(abs(dEyy0) < 1e-12);
            %dEyy0(indzeros) = 0;
            if useautocorr
                dsy = 1/2./sdy.*diag(dEyy0);
                dsy = dsy*sdy'+sdy*dsy';
                dEyy0corr = (dEyy0.*sy-dsy.*Eyy0)./(sy.*sy);
                dEyy0corr = dEyy0corr-diag(diag(dEyy0corr))+diag(diag(dEyy0));
                dMOMENTS(obs_nbr+1 : obs_nbr+obs_nbr*(obs_nbr+1)/2 , jp) = dyn_vech(dEyy0corr); %focus only on VAROBS variables
            else
                dMOMENTS(obs_nbr+1 : obs_nbr+obs_nbr*(obs_nbr+1)/2 , jp) = dyn_vech(dEyy0); %focus only on VAROBS variables
            end
            tmpEyyi = A*Ezz0*C(stationary_vars,:)' + B*Varinov*D(stationary_vars,:)';
            %we could distinguish between stderr and corr params, but this has no real speed effect as we multipliy with zeros
            dtmpEyyi = dA(:,:,jp)*Ezz0*C' + A*dEzz0*C' + A*Ezz0*dC(:,:,jp)' + dB(:,:,jp)*Varinov*D' + B*dVarinov(:,:,jp)*D' + B*Varinov*dD(:,:,jp)';
            Ai = eye(size(A,1)); %this is A^0
            dAi = zeros(size(A,1),size(A,1)); %this is d(A^0)
            for i = 1:nlags
                Eyyi = C(stationary_vars,:)*Ai*tmpEyyi;
                dEyyi = dC(:,:,jp)*Ai*tmpEyyi + C*dAi*tmpEyyi + C*Ai*dtmpEyyi;
                if useautocorr
                    dEyyi = (dEyyi.*sy-dsy.*Eyyi)./(sy.*sy);
                end
                dMOMENTS(obs_nbr + obs_nbr*(obs_nbr+1)/2 + (i-1)*obs_nbr^2 + 1 : obs_nbr + obs_nbr*(obs_nbr+1)/2 + i*obs_nbr^2, jp) = vec(dEyyi); %focus only on VAROBS variables
                dAi = dAi*A + Ai*dA(:,:,jp); %note that this is d(A^(i-1))
                Ai = Ai*A; %note that this is A^(i-1)
            end
        end
    end
else
    MOMENTS = [];
    dMOMENTS = [];
end

%% Compute dSPECTRUM
%Note that our state space system for computing the spectrum is the following:
% zhat = A*zhat(-1) + B*xi, where zhat = z - E(z)
% yhat = C*zhat(-1) + D*xi, where yhat = y - E(y)
if ~no_identification_spectrum
    %Some info on the spectral density: Dynare's state space system is given for states (z_{t} = A*z_{t-1} + B*u_{t}) and observables (y_{t} = C*z_{t-1} + D*u_{t})
        %The spectral density for y_{t} can be computed using different methods, which are numerically equivalent
        %See the following code example for the different ways;
%             freqs = 2*pi*(-(grid_nbr/2):1:(grid_nbr/2))'/grid_nbr; %divides the interval [-pi;pi] into ngrid+1 points
%             tpos  = exp( sqrt(-1)*freqs); %positive Fourier frequencies
%             tneg  = exp(-sqrt(-1)*freqs); %negative Fourier frequencies
%             IA = eye(size(A,1));
%             IE = eye(exo_nbr);
%             mathp_col1 = NaN(length(freqs),obs_nbr^2); mathp_col2 = mathp_col1; mathp_col3 = mathp_col1; mathp_col4 = mathp_col1;
%             for ig = 1:length(freqs)
%                 %method 1: as in UnivariateSpectralDensity.m
%                 f_omega  =(1/(2*pi))*( [(IA-A*tneg(ig))\B;IE]*Sigma_e*[B'/(IA-A'*tpos(ig)) IE]); % state variables
%                 g_omega1 = [C*tneg(ig) D]*f_omega*[C'*tpos(ig); D']; % selected variables
%                 %method 2: as in UnivariateSpectralDensity.m but simplified algebraically
%                 g_omega2 = (1/(2*pi))*(  C*((tpos(ig)*IA-A)\(B*Sigma_e*B'))*((tneg(ig)*IA-A')\(C'))  +  D*Sigma_e*B'*((tneg(ig)*IA-A')\(C'))  +   C* ((tpos(ig)*IA-A)\(B*Sigma_e*D'))  +  D*Sigma_e*D'  );
%                 %method 3: use transfer function note that ' is the complex conjugate transpose operator i.e. transpose(ffneg')==ffpos
%                 Transferfct = D+C*((tpos(ig)*IA-A)\B);
%                 g_omega3 = (1/(2*pi))*(Transferfct*Sigma_e*Transferfct');
%                 %method 4: kronecker products
%                 g_omega4 = (1/(2*pi))*( kron( D+C*((tneg(ig)^(-1)*IA-A)\B) , D+C*((tneg(ig)*IA-A)\B) )*Sigma_e(:));
%                 % store as matrix row
%                 mathp_col1(ig,:) = (g_omega1(:))'; mathp_col2(ig,:) = (g_omega2(:))'; mathp_col3(ig,:) = (g_omega3(:))'; mathp_col4(ig,:) = g_omega4;
%             end
%         disp([norm(mathp_col1 - mathp_col2); norm(mathp_col1 - mathp_col3); norm(mathp_col1 - mathp_col4); norm(mathp_col2 - mathp_col3); norm(mathp_col2 - mathp_col4); norm(mathp_col3 - mathp_col4);])
        %In the following we focus on method 3
    %Symmetry:
    %  Note that for the compuation of the G matrix we focus only on positive Fourier frequencies due to symmetry of the real part of the spectral density and, hence, the G matrix (which is real by construction).
    %  E.g. if grid_nbr=4, then we subdivide the intervall [-pi;pi] into [-3.1416;-1.5708;0;1.5708;3.1416], but focus only on [0;1.5708;3.1416] for the computation of the G matrix, 
    %  keeping in mind that the frequencies [1.5708;3.1416] need to be added twice, whereas the 0 frequency is only added once.    
    freqs = (0 : pi/(grid_nbr/2):pi); % we focus only on positive frequencies
    tpos  = exp( sqrt(-1)*freqs); %positive Fourier frequencies
    tneg  = exp(-sqrt(-1)*freqs); %negative Fourier frequencies
    IA = eye(size(A,1));
    if kronflag == -1
        %numerical derivative of spectral density
        dOmega_tmp = fjaco(str2func('identification_numerical_objective'), xparam1, 2, estim_params, M, oo, options, indpmodel, indpstderr, indpcorr, indvobs, useautocorr, nlags, grid_nbr); %[outputflag=2]
        kk = 0;
        for ig = 1:length(freqs)
            kk = kk+1;
            dOmega = dOmega_tmp(1 + (kk-1)*obs_nbr^2 : kk*obs_nbr^2,:);
            if ig == 1 % add zero frequency once
                dSPECTRUM = dOmega'*dOmega;
            else % due to symmetry to negative frequencies we can add positive frequencies twice
                dSPECTRUM = dSPECTRUM + 2*(dOmega'*dOmega);
            end
        end
    elseif kronflag == 1
        %use Kronecker products
        dA = reshape(dA,size(dA,1)*size(dA,2),size(dA,3));
        dB = reshape(dB,size(dB,1)*size(dB,2),size(dB,3));
        dC = reshape(dC,size(dC,1)*size(dC,2),size(dC,3));
        dD = reshape(dD,size(dD,1)*size(dD,2),size(dD,3));
        dVarinov = reshape(dVarinov,size(dVarinov,1)*size(dVarinov,2),size(dVarinov,3));
        K_obs_exo = commutation(obs_nbr,size(Varinov,1));
        for ig=1:length(freqs)
            z = tneg(ig);
            zIminusA =  (z*IA - A);
            zIminusAinv = zIminusA\IA;
            Transferfct = D + C*zIminusAinv*B; % Transfer function
            dzIminusA = -dA;
            dzIminusAinv = kron(-(transpose(zIminusA)\IA),zIminusAinv)*dzIminusA; %this takes long
            dTransferfct = dD + DerivABCD(C,dC,zIminusAinv,dzIminusAinv,B,dB); %also long
            dTransferfct_conjt = K_obs_exo*conj(dTransferfct);
            dOmega = (1/(2*pi))*DerivABCD(Transferfct,dTransferfct,Varinov,dVarinov,Transferfct',dTransferfct_conjt); %also long
            if ig == 1 % add zero frequency once
                dSPECTRUM = dOmega'*dOmega;
            else  % due to symmetry to negative frequencies we can add positive frequencies twice
                dSPECTRUM = dSPECTRUM + 2*(dOmega'*dOmega);
            end
        end
        %put back into tensor notation
        dA = reshape(dA,size(A,1),size(A,2),totparam_nbr);
        dB = reshape(dB,size(B,1),size(B,2),totparam_nbr);
        dC = reshape(dC,size(C,1),size(C,2),totparam_nbr);
        dD = reshape(dD,size(D,1),size(D,2),totparam_nbr);
        dVarinov = reshape(dVarinov,size(Varinov,1),size(Varinov,2),totparam_nbr);
    elseif (kronflag==0) || (kronflag==-2)
        for ig = 1:length(freqs)
            IzminusA =  tpos(ig)*IA - A;
            invIzminusA = IzminusA\eye(size(A,1));
            Transferfct = D + C*invIzminusA*B;
            dOmega = zeros(obs_nbr^2,totparam_nbr);
            for j = 1:totparam_nbr
                if j <= stderrparam_nbr+corrparam_nbr %stderr and corr parameters: only dSig is nonzero
                    dOmega_tmp = Transferfct*dVarinov(:,:,j)*Transferfct';
                else %model parameters
                    dinvIzminusA = -invIzminusA*(-dA(:,:,j))*invIzminusA;
                    dTransferfct = dD(:,:,j) + dC(:,:,j)*invIzminusA*B + C*dinvIzminusA*B + C*invIzminusA*dB(:,:,j);
                    dOmega_tmp = dTransferfct*Varinov*Transferfct' + Transferfct*dVarinov(:,:,j)*Transferfct' + Transferfct*Varinov*dTransferfct';
                end
                dOmega(:,j) = (1/(2*pi))*dOmega_tmp(:);
            end
            if ig == 1 % add zero frequency once
                dSPECTRUM = dOmega'*dOmega;
            else % due to symmetry to negative frequencies we can add positive frequencies twice
                dSPECTRUM = dSPECTRUM + 2*(dOmega'*dOmega);
            end
        end
    end
    % Normalize Matrix and add steady state Jacobian, note that G is real and symmetric by construction
    dSPECTRUM = 2*pi*dSPECTRUM./(2*length(freqs)-1) + dMEAN'*dMEAN;
    dSPECTRUM = real(dSPECTRUM);
else
    dSPECTRUM = [];
end

%% Compute dMINIMAL
if ~no_identification_minimal
    if obs_nbr < exo_nbr
        % Check whether criteria can be used
        warning_KomunjerNg = 'WARNING: Komunjer and Ng (2011) failed:\n';
        warning_KomunjerNg = [warning_KomunjerNg '       There are more shocks and measurement errors than observables, this is not implemented (yet).\n'];
        warning_KomunjerNg = [warning_KomunjerNg '       Skip identification analysis based on minimal state space system.\n'];
        fprintf(warning_KomunjerNg);        
        dMINIMAL = [];        
    else
        % Derive and check minimal state vector of first-order
        [CheckCO,minnx,minA,minB,minC,minD,dminA,dminB,dminC,dminD] = get_minimal_state_representation(oo.dr.ghx(oo.dr.pruned.indx,:),...            %A
                                                                                                       oo.dr.ghu(oo.dr.pruned.indx,:),...            %B
                                                                                                       oo.dr.ghx(oo.dr.pruned.indy,:),...             %C
                                                                                                       oo.dr.ghu(oo.dr.pruned.indy,:),...             %D
                                                                                                       oo.dr.derivs.dghx(oo.dr.pruned.indx,:,:),...  %dA
                                                                                                       oo.dr.derivs.dghu(oo.dr.pruned.indx,:,:),...  %dB
                                                                                                       oo.dr.derivs.dghx(oo.dr.pruned.indy,:,:),...   %dC
                                                                                                       oo.dr.derivs.dghu(oo.dr.pruned.indy,:,:));     %dD
        if CheckCO == 0
            warning_KomunjerNg = 'WARNING: Komunjer and Ng (2011) failed:\n';
            warning_KomunjerNg = [warning_KomunjerNg '         Conditions for minimality are not fullfilled:\n'];
            warning_KomunjerNg = [warning_KomunjerNg '         Skip identification analysis based on minimal state space system.\n'];
            fprintf(warning_KomunjerNg); %use sprintf to have line breaks            
            dMINIMAL = [];
        else
            %reshape into Magnus-Neudecker Jacobians, i.e. dvec(X)/dp
            dminA = reshape(dminA,size(dminA,1)*size(dminA,2),size(dminA,3));
            dminB = reshape(dminB,size(dminB,1)*size(dminB,2),size(dminB,3));
            dminC = reshape(dminC,size(dminC,1)*size(dminC,2),size(dminC,3));
            dminD = reshape(dminD,size(dminD,1)*size(dminD,2),size(dminD,3));
            dvechSig = reshape(oo.dr.derivs.dSigma_e,exo_nbr*exo_nbr,totparam_nbr);
            indvechSig= find(tril(ones(exo_nbr,exo_nbr)));
            dvechSig = dvechSig(indvechSig,:);
            Inx = eye(minnx);
            Inu = eye(exo_nbr);
            [~,Enu] = duplication(exo_nbr);
            KomunjerNg_DL = [dminA; dminB; dminC; dminD; dvechSig];
            KomunjerNg_DT = [kron(transpose(minA),Inx) - kron(Inx,minA);
                             kron(transpose(minB),Inx);
                             -1*kron(Inx,minC);
                             zeros(obs_nbr*exo_nbr,minnx^2);
                             zeros(exo_nbr*(exo_nbr+1)/2,minnx^2)];
            KomunjerNg_DU = [zeros(minnx^2,exo_nbr^2);
                             kron(Inu,minB);
                             zeros(obs_nbr*minnx,exo_nbr^2);
                             kron(Inu,minD);
                             -2*Enu*kron(Sigma_e,Inu)];
            dMINIMAL = full([KomunjerNg_DL KomunjerNg_DT KomunjerNg_DU]);
            %add Jacobian of steady state (here we also allow for higher-order perturbation, i.e. only the mean provides additional restrictions
            dMINIMAL =  [dMEAN zeros(obs_nbr,minnx^2+exo_nbr^2); dMINIMAL];
        end
    end
else
    dMINIMAL = [];
end


function [dX] = DerivABCD(X1,dX1,X2,dX2,X3,dX3,X4,dX4)
% function [dX] = DerivABCD(X1,dX1,X2,dX2,X3,dX3,X4,dX4)
% -------------------------------------------------------------------------
% Derivative of X(p)=X1(p)*X2(p)*X3(p)*X4(p) w.r.t to p
% See Magnus and Neudecker (1999), p. 175
% -------------------------------------------------------------------------
% Inputs: Matrices X1,X2,X3,X4, and the corresponding derivatives w.r.t p.
% Output: Derivative of product of X1*X2*X3*X4 w.r.t. p
% =========================================================================
nparam = size(dX1,2);
%   If one or more matrices are left out, they are set to zero
if nargin == 4
    X3=speye(size(X2,2)); dX3=spalloc(numel(X3),nparam,0);
    X4=speye(size(X3,2)); dX4=spalloc(numel(X4),nparam,0);
elseif nargin == 6
    X4=speye(size(X3,2)); dX4=spalloc(numel(X4),nparam,0);
end

dX = kron(transpose(X4)*transpose(X3)*transpose(X2),speye(size(X1,1)))*dX1...
   + kron(transpose(X4)*transpose(X3),X1)*dX2...
   + kron(transpose(X4),X1*X2)*dX3...
   + kron(speye(size(X4,2)),X1*X2*X3)*dX4;
end %DerivABCD end

end%main function end