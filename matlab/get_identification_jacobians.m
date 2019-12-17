function [J, G, D, dTAU, dLRE, dA, dOm, dYss, MOMENTS] = get_identification_jacobians(A, B, estim_params, M, oo, options, options_ident, indpmodel, indpstderr, indpcorr, indvobs)
% [J, G, D, dTAU, dLRE, dA, dOm, dYss, MOMENTS] = get_identification_jacobians(A, B, estim_params, M, oo, options, options_ident, indpmodel, indpstderr, indpcorr, indvobs)
% previously getJJ.m
% -------------------------------------------------------------------------
% Sets up the Jacobians needed for identification analysis based on the 
% Iskrev's J, Qu and Tkachenko's G and Komunjer and Ng's D matrices as well
% as on the reduced-form model (dTAU) and the dynamic model (dLRE)
% =========================================================================
% INPUTS
%   A:              [endo_nbr by endo_nbr] in DR order
%                     Transition matrix from Kalman filter for all endogenous declared variables, 
%   B:              [endo_nbr by exo_nbr] in DR order
%                     Transition matrix from Kalman filter mapping shocks today to endogenous variables today
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
%   J:              [obs_nbr+obs_nbr*(obs_nbr+1)/2+nlags*obs_nbr^2 by totparam_nbr] in DR Order
%                     Jacobian of 1st and 2nd order moments of observables wrt, 
%                     all parameters, i.e. dgam (see below for definition of gam).
%                     Corresponds to Iskrev (2010)'s J matrix.
%   G:              [totparam_nbr by totparam_nbr] in DR Order
%                     Sum of (1) Gram Matrix of Jacobian of spectral density
%                     wrt to all parameters plus (2) Gram Matrix of Jacobian
%                     of steady state wrt to all parameters.
%                     Corresponds to Qu and Tkachenko (2012)'s G matrix.
%   D:              [obs_nbr+minstate_nbr^2+minstate_nbr*exo_nbr+obs_nbr*minstate_nbr+obs_nbr*exo_nbr+exo_nbr*(exo_nbr+1)/2  by  totparam_nbr+minstate_nbr^2+exo_nbr^2] in DR order
%                     Jacobian of minimal System Matrices and unique Transformation 
%                     of steady state and spectral density matrix.
%                     Corresponds to Komunjer and Ng (2011)'s Delta matrix.
%   dTAU:           [(endo_nbr+endo_nbr^2+endo_nbr*(endo_nbr+1)/2) by totparam_nbr] in DR order
%                     Jacobian of linearized reduced form state space model, given Yss [steady state],
%                     A [transition matrix], B [matrix of shocks], Sigma [covariance of shocks]
%                     tau = [ys; vec(A); dyn_vech(B*Sigma*B')] with respect to all parameters.
%   dLRE:           [endo_nbr+endo_nbr*(dynamicvar_nbr+exo_nbr) by modparam_nbr] in DR order
%                     Jacobian of steady state and linear rational expectation matrices
%                     (i.e. Jacobian of dynamic model) with respect to estimated model parameters only (indpmodel)
%   dA:             [endo_nbr by endo_nbr by totparam_nbr] in DR order
%                     Jacobian (wrt to all parameters) of transition matrix A
%   dOm:            [endo_nbr by endo_nbr by totparam_nbr] in DR order
%                     Jacobian (wrt to all paramters) of Om = (B*Sigma_e*B')
%   dYss            [endo_nbr by modparam_nbr] in DR order
%                     Jacobian (wrt model parameters only) of steady state
%   MOMENTS:        [obs_nbr+obs_nbr*(obs_nbr+1)/2+nlags*obs_nbr^2 by 1]
%                     vector of theoretical moments of observed (VAROBS)
%                     variables. Note that J is the Jacobian of MOMENTS.
%                     MOMENTS = [ys(indvobs); dyn_vech(GAM{1}); vec(GAM{j+1})]; for j=1:ar and where
%                     GAM is the first output of th_autocovariances
% -------------------------------------------------------------------------
% This function is called by 
%   * identification_analysis.m
% -------------------------------------------------------------------------
% This function calls
%   * commutation
%   * get_minimal_state_representation
%   * duplication
%   * dyn_vech
%   * get_perturbation_params_deriv (previously getH)
%   * get_all_parameters
%   * fjaco
%   * lyapunov_symm
%   * th_autocovariances
%   * identification_numerical_objective (previously thet2tau)
%   * vec
% =========================================================================
% Copyright (C) 2010-2019 Dynare Team
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

%get options
nlags       = options_ident.ar;
useautocorr = options_ident.useautocorr;
grid_nbr    = options_ident.grid_nbr;
kronflag    = options_ident.analytic_derivation_mode;
no_identification_moments   = options_ident.no_identification_moments;
no_identification_minimal   = options_ident.no_identification_minimal;
no_identification_spectrum  = options_ident.no_identification_spectrum;
params0     = M.params;                 %values at which to evaluate dynamic, static and param_derivs files
Sigma_e0    = M.Sigma_e;                %covariance matrix of exogenous shocks
Corr_e0     = M.Correlation_matrix;     %correlation matrix of exogenous shocks
stderr_e0   = sqrt(diag(Sigma_e0));     %standard errors of exogenous shocks
para0       = get_all_parameters(estim_params, M);  %get all selected parameters in estimated_params block, stderr and corr come first, then model parameters
if isempty(para0)
    %if there is no estimated_params block, consider all stderr and all model parameters, but no corr parameters
    para0 = [stderr_e0', params0'];
end
%get numbers/lengths of vectors
modparam_nbr    = length(indpmodel);
stderrparam_nbr = length(indpstderr);
corrparam_nbr   = size(indpcorr,1);
totparam_nbr    = stderrparam_nbr + corrparam_nbr + modparam_nbr;
obs_nbr         = length(indvobs);
exo_nbr         = M.exo_nbr;
endo_nbr        = M.endo_nbr;

%% Construct dTAU, dLRE, dA, dOm, dYss
DERIVS = get_perturbation_params_derivs(M, options, estim_params, oo, indpmodel, indpstderr, indpcorr, 0);%d2flag=0
dA = zeros(M.endo_nbr, M.endo_nbr, size(DERIVS.dghx,3));
dA(:,M.nstatic+(1:M.nspred),:) = DERIVS.dghx;
dB = DERIVS.dghu;
dSig = DERIVS.dSigma_e;
dOm = DERIVS.dOm;
dYss = DERIVS.dYss;

% Collect terms for derivative of tau=[Yss; vec(A); vech(Om)]
dTAU = zeros(endo_nbr*endo_nbr+endo_nbr*(endo_nbr+1)/2, totparam_nbr);
for j=1:totparam_nbr
    dTAU(:,j) = [vec(dA(:,:,j)); dyn_vech(dOm(:,:,j))];
end
dTAU = [ [zeros(endo_nbr, stderrparam_nbr+corrparam_nbr) dYss]; dTAU ]; % add steady state
dLRE = [dYss; reshape(DERIVS.dg1,size(DERIVS.dg1,1)*size(DERIVS.dg1,2),size(DERIVS.dg1,3)) ]; %reshape dg1 and add steady state
    
%State Space Matrices for VAROBS variables
C  = A(indvobs,:);
dC = dA(indvobs,:,:);
D  = B(indvobs,:);
dD = dB(indvobs,:,:);
    
%% Iskrev (2010)
if ~no_identification_moments
    if kronflag == -1
        %numerical derivative of autocovariogram [outputflag=1]
        J = fjaco(str2func('identification_numerical_objective'), para0, 1, estim_params, M, oo, options, indpmodel, indpstderr, indpcorr, indvobs, useautocorr, nlags, grid_nbr);
        M.params = params0;              %make sure values are set back
        M.Sigma_e = Sigma_e0;            %make sure values are set back
        M.Correlation_matrix = Corr_e0 ; %make sure values are set back
        J = [[zeros(obs_nbr,stderrparam_nbr+corrparam_nbr) dYss(indvobs,:)]; J]; %add Jacobian of steady state of VAROBS variables
    else
        J = zeros(obs_nbr + obs_nbr*(obs_nbr+1)/2 + nlags*obs_nbr^2 , totparam_nbr);
        J(1:obs_nbr,stderrparam_nbr+corrparam_nbr+1 : totparam_nbr) = dYss(indvobs,:); %add Jacobian of steady state of VAROBS variables
        % Denote Ezz0 = E_t(z_t * z_t'), then the following Lyapunov equation defines the autocovariagram: Ezz0 -A*Ezz*A' = B*Sig_e*B' = Om        
        Gamma_y =  lyapunov_symm(A, B*Sigma_e0*B', options.lyapunov_fixed_point_tol, options.qz_criterium, options.lyapunov_complex_threshold, 1, options.debug);
        %here method=1 is used, whereas all other calls of lyapunov_symm use method=2. The reason is that T and U are persistent, and input matrix A is the same, so using option 2 for all the rest of iterations spares a lot of computing time while not repeating Schur every time
        indzeros = find(abs(Gamma_y) < 1e-12); %find values that are numerical zero
        Gamma_y(indzeros) = 0;
        %   if useautocorr,
        sdy = sqrt(diag(Gamma_y)); %theoretical standard deviation
        sy = sdy*sdy';             %cross products of standard deviations
        %   end
        for j = 1:(stderrparam_nbr+corrparam_nbr)
            %Jacobian of Ezz0 wrt exogenous paramters: dEzz0(:,:,j)-A*dEzz0(:,:,j)*A'=dOm(:,:,j), because dA is zero by construction for stderr and corr parameters
            dEzz0 =  lyapunov_symm(A,dOm(:,:,j),options.lyapunov_fixed_point_tol,options.qz_criterium,options.lyapunov_complex_threshold,2,options.debug);
            %here method=2 is used to spare a lot of computing time while not repeating Schur every time
            indzeros = find(abs(dEzz0) < 1e-12);
            dEzz0(indzeros) = 0;
            if useautocorr
                dsy = 1/2./sdy.*diag(dEzz0);
                dsy = dsy*sdy'+sdy*dsy';
                dEzz0corr = (dEzz0.*sy-dsy.*Gamma_y)./(sy.*sy);
                dEzz0corr = dEzz0corr-diag(diag(dEzz0corr))+diag(diag(dEzz0));
                J(obs_nbr+1 : obs_nbr+obs_nbr*(obs_nbr+1)/2 , j) = dyn_vech(dEzz0corr(indvobs,indvobs)); %focus only on VAROBS variables
            else
                J(obs_nbr+1 : obs_nbr+obs_nbr*(obs_nbr+1)/2 , j) = dyn_vech(dEzz0(indvobs,indvobs)); %focus only on VAROBS variables
            end
            %Jacobian of Ezzi = E_t(z_t * z_{t-i}'): dEzzi(:,:,j) = A^i*dEzz0(:,:,j) wrt stderr and corr parameters, because dA is zero by construction for stderr and corr parameters
            for i = 1:nlags
                dEzzi = A^i*dEzz0;
                if useautocorr
                    dEzzi = (dEzzi.*sy-dsy.*(A^i*Gamma_y))./(sy.*sy);
                end
                J(obs_nbr + obs_nbr*(obs_nbr+1)/2 + (i-1)*obs_nbr^2 + 1 : obs_nbr + obs_nbr*(obs_nbr+1)/2 + i*obs_nbr^2, j) = vec(dEzzi(indvobs,indvobs)); %focus only on VAROBS variables
            end
        end
        for j=1:modparam_nbr
            %Jacobian of Ezz0 wrt model parameters: dEzz0(:,:,j) - A*dEzz0(:,:,j)*A' = dOm(:,:,j) + dA(:,:,j)*Ezz*A'+ A*Ezz*dA(:,:,j)'
            dEzz0 =  lyapunov_symm(A,dA(:,:,j+stderrparam_nbr+corrparam_nbr)*Gamma_y*A'+A*Gamma_y*dA(:,:,j+stderrparam_nbr+corrparam_nbr)'+dOm(:,:,j+stderrparam_nbr+corrparam_nbr),options.lyapunov_fixed_point_tol,options.qz_criterium,options.lyapunov_complex_threshold,2,options.debug);
            %here method=2 is used to spare a lot of computing time while not repeating Schur every time
            indzeros = find(abs(dEzz0) < 1e-12);
            dEzz0(indzeros) = 0;
            if useautocorr
                dsy = 1/2./sdy.*diag(dEzz0);
                dsy = dsy*sdy'+sdy*dsy';
                dEzz0corr = (dEzz0.*sy-dsy.*Gamma_y)./(sy.*sy);
                dEzz0corr = dEzz0corr-diag(diag(dEzz0corr))+diag(diag(dEzz0));
                J(obs_nbr+1 : obs_nbr+obs_nbr*(obs_nbr+1)/2 , stderrparam_nbr+corrparam_nbr+j) = dyn_vech(dEzz0corr(indvobs,indvobs)); %focus only on VAROBS variables
            else
                J(obs_nbr+1 : obs_nbr+obs_nbr*(obs_nbr+1)/2 , stderrparam_nbr+corrparam_nbr+j) = dyn_vech(dEzz0(indvobs,indvobs)); %focus only on VAROBS variables
            end
            %Jacobian of Ezzi = E_t(z_t * z_{t-i}'): dEzzi(:,:,j) = A^i*dEzz0(:,:,j) + d(A^i)*dEzz0(:,:,j) wrt model parameters        
            for i = 1:nlags            
                dEzzi = A^i*dEzz0;
                for ii=1:i
                    dEzzi = dEzzi + A^(ii-1)*dA(:,:,j+stderrparam_nbr+corrparam_nbr)*A^(i-ii)*Gamma_y;
                end
                if useautocorr
                    dEzzi = (dEzzi.*sy-dsy.*(A^i*Gamma_y))./(sy.*sy);
                end
                J(obs_nbr + obs_nbr*(obs_nbr+1)/2 + (i-1)*obs_nbr^2 + 1 : obs_nbr + obs_nbr*(obs_nbr+1)/2 + i*obs_nbr^2, j+stderrparam_nbr+corrparam_nbr) = vec(dEzzi(indvobs,indvobs)); %focus only on VAROBS variables
            end
        end
    end
else
    J = [];
end
    
%% Qu and Tkachenko (2012)
if ~no_identification_spectrum
    %Some info on the spectral density: Dynare's state space system is given for states (z_{t} = A*z_{t-1} + B*u_{t}) and observables (y_{t} = C*z_{t-1} + D*u_{t})
        %The spectral density for y_{t} can be computed using different methods, which are numerically equivalent
        %See the following code example for the different ways;
%             freqs = 2*pi*(-(grid_nbr/2):1:(grid_nbr/2))'/grid_nbr; %divides the interval [-pi;pi] into ngrid+1 points
%             tpos  = exp( sqrt(-1)*freqs); %positive Fourier frequencies
%             tneg  = exp(-sqrt(-1)*freqs); %negative Fourier frequencies
%             IA = eye(size(A,1));
%             IE = eye(M.exo_nbr);
%             mathp_col1 = NaN(length(freqs),obs_nbr^2); mathp_col2 = mathp_col1; mathp_col3 = mathp_col1; mathp_col4 = mathp_col1;
%             for ig = 1:length(freqs)
%                 %method 1: as in UnivariateSpectralDensity.m
%                 f_omega  =(1/(2*pi))*( [(IA-A*tneg(ig))\B;IE]*M.Sigma_e*[B'/(IA-A'*tpos(ig)) IE]); % state variables
%                 g_omega1 = [C*tneg(ig) D]*f_omega*[C'*tpos(ig); D']; % selected variables
%                 %method 2: as in UnivariateSpectralDensity.m but simplified algebraically
%                 g_omega2 = (1/(2*pi))*(  C*((tpos(ig)*IA-A)\(B*M.Sigma_e*B'))*((tneg(ig)*IA-A')\(C'))  +  D*M.Sigma_e*B'*((tneg(ig)*IA-A')\(C'))  +   C* ((tpos(ig)*IA-A)\(B*M.Sigma_e*D'))  +  D*M.Sigma_e*D'  );
%                 %method 3: use transfer function note that ' is the complex conjugate transpose operator i.e. transpose(ffneg')==ffpos
%                 Transferfct = D+C*((tpos(ig)*IA-A)\B);
%                 g_omega3 = (1/(2*pi))*(Transferfct*M.Sigma_e*Transferfct');
%                 %method 4: kronecker products
%                 g_omega4 = (1/(2*pi))*( kron( D+C*((tneg(ig)^(-1)*IA-A)\B) , D+C*((tneg(ig)*IA-A)\B) )*M.Sigma_e(:));
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
        %numerical derivative of spectral density [outputflag=2]
        dOmega_tmp = fjaco(str2func('identification_numerical_objective'), para0, 2, estim_params, M, oo, options, indpmodel, indpstderr, indpcorr, indvobs, useautocorr, nlags, grid_nbr);
        M.params = params0;              %make sure values are set back
        M.Sigma_e = Sigma_e0;            %make sure values are set back
        M.Correlation_matrix = Corr_e0 ; %make sure values are set back
        kk = 0;
        for ig = 1:length(freqs)
            kk = kk+1;
            dOmega = dOmega_tmp(1 + (kk-1)*obs_nbr^2 : kk*obs_nbr^2,:);
            if ig == 1 % add zero frequency once
                G = dOmega'*dOmega;
            else % due to symmetry to negative frequencies we can add positive frequencies twice
                G = G + 2*(dOmega'*dOmega);
            end
        end
    elseif kronflag == 1        
        %use Kronecker products
        dA = reshape(dA,size(dA,1)*size(dA,2),size(dA,3));
        dB = reshape(dB,size(dB,1)*size(dB,2),size(dB,3));
        dC = reshape(dC,size(dC,1)*size(dC,2),size(dC,3));
        dD = reshape(dD,size(dD,1)*size(dD,2),size(dD,3));
        dSig = reshape(dSig,size(dSig,1)*size(dSig,2),size(dSig,3));
        K_obs_exo = commutation(obs_nbr,exo_nbr);
        for ig=1:length(freqs)
            z = tneg(ig);
            zIminusA =  (z*IA - A);
            zIminusAinv = zIminusA\IA;
            Transferfct = D + C*zIminusAinv*B; % Transfer function
            dzIminusA = -dA;
            dzIminusAinv = kron(-(transpose(zIminusA)\IA),zIminusAinv)*dzIminusA;
            dTransferfct = dD + DerivABCD(C,dC,zIminusAinv,dzIminusAinv,B,dB);
            dTransferfct_conjt = K_obs_exo*conj(dTransferfct);
            dOmega = (1/(2*pi))*DerivABCD(Transferfct,dTransferfct,Sigma_e0,dSig,Transferfct',dTransferfct_conjt);
            if ig == 1 % add zero frequency once
                G = dOmega'*dOmega;
            else  % due to symmetry to negative frequencies we can add positive frequencies twice
                G = G + 2*(dOmega'*dOmega);
            end
        end
        %put back into tensor notation
        dA = reshape(dA,endo_nbr,endo_nbr,totparam_nbr);
        dB = reshape(dB,endo_nbr,exo_nbr,totparam_nbr);
        dC = reshape(dC,obs_nbr,endo_nbr,totparam_nbr);
        dD = reshape(dD,obs_nbr,exo_nbr,totparam_nbr);
        dSig = reshape(dSig,exo_nbr,exo_nbr,totparam_nbr);
    elseif (kronflag==0) || (kronflag==-2)
        for ig = 1:length(freqs)
            IzminusA =  tpos(ig)*IA - A;            
            invIzminusA = IzminusA\eye(endo_nbr);
            Transferfct = D + C*invIzminusA*B;
            dOmega = zeros(obs_nbr^2,totparam_nbr);
            for j = 1:totparam_nbr
                if j <= stderrparam_nbr+corrparam_nbr %stderr and corr parameters: only dSig is nonzero
                    dOmega_tmp = Transferfct*dSig(:,:,j)*Transferfct';
                else %model parameters
                    dinvIzminusA = -invIzminusA*(-dA(:,:,j))*invIzminusA;
                    dTransferfct = dD(:,:,j) + dC(:,:,j)*invIzminusA*B + C*dinvIzminusA*B + C*invIzminusA*dB(:,:,j);
                    dOmega_tmp = dTransferfct*M.Sigma_e*Transferfct' + Transferfct*dSig(:,:,j)*Transferfct' + Transferfct*M.Sigma_e*dTransferfct';
                end
                dOmega(:,j) = (1/(2*pi))*dOmega_tmp(:);
            end
            if ig == 1 % add zero frequency once
                G = dOmega'*dOmega;
            else % due to symmetry to negative frequencies we can add positive frequencies twice
                G = G + 2*(dOmega'*dOmega);
            end
        end        
    end    
    % Normalize Matrix and add steady state Jacobian, note that G is real and symmetric by construction
    G = 2*pi*G./(2*length(freqs)-1) + [zeros(obs_nbr,stderrparam_nbr+corrparam_nbr) dYss(indvobs,:)]'*[zeros(obs_nbr,stderrparam_nbr+corrparam_nbr) dYss(indvobs,:)];
    G = real(G);
else
    G = [];
end

%% Komunjer and Ng (2012)
if ~no_identification_minimal
    if obs_nbr < exo_nbr
        % Check whether criteria can be used
        warning_KomunjerNg = 'WARNING: Komunjer and Ng (2011) failed:\n';
        warning_KomunjerNg = [warning_KomunjerNg '       There are more shocks and measurement errors than observables, this is not implemented (yet).\n'];
        warning_KomunjerNg = [warning_KomunjerNg '       Skip identification analysis based on minimal state space system.\n'];
        fprintf(warning_KomunjerNg);        
        D = [];        
    else
        % Derive and check minimal state        
        if isfield(oo.dr,'state_var')            
            state_var = oo.dr.state_var; %state variables in declaration order
        else
            % DR-order: static variables first, then purely backward variables, then mixed variables, finally purely forward variables.
            % Inside each category, variables are arranged according to the declaration order.
            % state variables are the purely backward variables and the mixed variables
            state_var = transpose(oo.dr.order_var( (M.nstatic+1):(M.nstatic+M.npred+M.nboth) ) ); %state variables in declaration order            
        end
        state_var_DR = oo.dr.inv_order_var(state_var); %state vector in DR order
        minA = A(state_var_DR,state_var_DR); dminA = dA(state_var_DR,state_var_DR,:);
        minB = B(state_var_DR,:);            dminB = dB(state_var_DR,:,:);
        minC = C(:,state_var_DR);            dminC = dC(:,state_var_DR,:);
        minD = D(:,:);                       dminD = dD(:,:,:);
        [CheckCO,minnx,minA,minB,minC,minD,dminA,dminB,dminC,dminD] = get_minimal_state_representation(minA,minB,minC,minD,dminA,dminB,dminC,dminD);        
        if CheckCO == 0
            warning_KomunjerNg = 'WARNING: Komunjer and Ng (2011) failed:\n';
            warning_KomunjerNg = [warning_KomunjerNg '         Conditions for minimality are not fullfilled:\n'];
            warning_KomunjerNg = [warning_KomunjerNg '         Skip identification analysis based on minimal state space system.\n'];
            fprintf(warning_KomunjerNg); %use sprintf to have line breaks            
            D = [];
        else
            %reshape into Magnus-Neudecker Jacobians, i.e. dvec(X)/dp
            dminA = reshape(dminA,size(dminA,1)*size(dminA,2),size(dminA,3));
            dminB = reshape(dminB,size(dminB,1)*size(dminB,2),size(dminB,3));
            dminC = reshape(dminC,size(dminC,1)*size(dminC,2),size(dminC,3));
            dminD = reshape(dminD,size(dminD,1)*size(dminD,2),size(dminD,3));
            dvechSig = reshape(dSig,size(dSig,1)*size(dSig,2),size(dSig,3));
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
                             -2*Enu*kron(M.Sigma_e,Inu)];
            D = full([KomunjerNg_DL KomunjerNg_DT KomunjerNg_DU]);
            %add Jacobian of steady state
            D =  [[zeros(obs_nbr,stderrparam_nbr+corrparam_nbr) dYss(indvobs,:)], zeros(obs_nbr,minnx^2+exo_nbr^2); D];
        end
    end
else
    D = [];
end

if nargout > 8
    options.ar = nlags;
    nodecomposition = 1;
    [Gamma_y,~] = th_autocovariances(oo.dr,oo.dr.order_var(indvobs),M,options,nodecomposition); %focus only on observed variables
    sdy=sqrt(diag(Gamma_y{1})); %theoretical standard deviation
    sy=sdy*sdy';                %cross products of standard deviations
    if useautocorr
        sy=sy-diag(diag(sy))+eye(obs_nbr);
        Gamma_y{1}=Gamma_y{1}./sy;
    else
        for j=1:nlags
            Gamma_y{j+1}=Gamma_y{j+1}.*sy;
        end
    end
    %Ezz is vector of theoretical moments of VAROBS variables
    MOMENTS = dyn_vech(Gamma_y{1});
    for j=1:nlags
        MOMENTS = [MOMENTS; vec(Gamma_y{j+1})];
    end
    MOMENTS = [oo.dr.ys(oo.dr.order_var(indvobs)); MOMENTS];    
end


function [dX] = DerivABCD(A,dA,B,dB,C,dC,D,dD)
% function [dX] = DerivABCD(A,dA,B,dB,C,dC,D,dD)
% -------------------------------------------------------------------------
% Derivative of X(p)=A(p)*B(p)*C(p)*D(p) w.r.t to p
% See Magnus and Neudecker (1999), p. 175
% -------------------------------------------------------------------------
% Inputs: Matrices A,B,C,D, and the corresponding derivatives w.r.t p.
% Output: Derivative of product of A*B*C*D w.r.t. p
% =========================================================================
nparam = size(dA,2);
%   If one or more matrices are left out, they are set to zero
if nargin == 4
    C=speye(size(B,2)); dC=spalloc(numel(C),nparam,0);
    D=speye(size(C,2)); dD=spalloc(numel(D),nparam,0);
elseif nargin == 6
    D=speye(size(C,2)); dD=spalloc(numel(D),nparam,0);
end

dX1 = kron(transpose(D)*transpose(C)*transpose(B),speye(size(A,1)))*dA;
dX2 = kron(transpose(D)*transpose(C),A)*dB;
dX3 = kron(transpose(D),A*B)*dC;
dX4 = kron(speye(size(D,2)),A*B*C)*dD;
dX= dX1+dX2+dX3+dX4;
end %DerivABCD end

end%main function end


