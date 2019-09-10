function [ide_moments, ide_spectrum, ide_minimal, ide_hess, ide_reducedform, ide_lre, derivatives_info, info, options_ident] = identification_analysis(params, indpmodel, indpstderr, indpcorr, options_ident, dataset_info, prior_exist, init)
% [ide_moments, ide_spectrum, ide_minimal, ide_hess, ide_reducedform, ide_lre, derivatives_info, info, options_ident] = identification_analysis(params, indpmodel, indpstderr, indpcorr, options_ident, dataset_info, prior_exist, init)
% -------------------------------------------------------------------------
% This function wraps all identification analysis, i.e. it
% (1) wraps functions for the theoretical identification analysis based on moments (Iskrev, 2010),
% spectrum (Qu and Tkachenko, 2012), minimal system (Komunjer and Ng, 2011),  information matrix,
% reduced-form solution and linear rational expectation model (Ratto and Iskrev, 2011).
% (2) computes the identification strength based on moments (Ratto and Iskrev, 2011)
% (3) checks which parameters are involved.
% =========================================================================
% INPUTS
%    * params             [mc_sample_nbr by totparam_nbr]
%                         parameter values for identification checks
%    * indpmodel          [modparam_nbr by 1]
%                         index of model parameters for which identification is checked
%    * indpstderr         [stderrparam_nbr by 1]
%                         index of stderr parameters for which identification is checked
%    * indpcorr           [corrparam_nbr by 1]
%                         index of corr parmeters for which identification is checked
%    * options_ident      [structure]
%                         identification options
%    * dataset_info       [structure]
%                         various information about the dataset (descriptive statistics and missing observations) for Kalman Filter
%    * prior_exist        [integer]
%                         1: prior exists. Identification is checked for stderr, corr and model parameters as declared in estimated_params block
%                         0: prior does not exist. Identification is checked for all stderr and model parameters, but no corr parameters
%    * init               [integer]
%                         flag for initialization of persistent vars. This is needed if one may want to re-initialize persistent
%                         variables even if they are not empty, e.g. by making more calls to identification in the same dynare mod file
% -------------------------------------------------------------------------
% OUTPUTS
%    * ide_moments        [structure]
%                         identification results using theoretical moments (Iskrev, 2010)
%    * ide_spectrum       [structure]
%                         identification results using spectrum (Qu and Tkachenko, 2012)
%    * ide_minimal        [structure]
%                         identification results using minimal system (Komunjer and Ng, 2011)
%    * ide_hess           [structure]
%                         identification results using asymptotic Hessian (Ratto and Iskrev, 2011)
%    * ide_reducedform    [structure]
%                         identification results using reduced form solution (Ratto and Iskrev, 2011)
%    * ide_lre            [structure]
%                         identification results using linear rational expectations model,
%                         i.e. steady state and Jacobian of dynamic model, (Ratto and Iskrev, 2011)
%    * derivatives_info   [structure]
%                         info about derivatives, used in dsge_likelihood.m
%    * info               [integer]
%                         output from dynare_resolve
%    * options_ident      [structure]
%                         updated identification options
% -------------------------------------------------------------------------
% This function is called by
%   * dynare_identification.m
% -------------------------------------------------------------------------
% This function calls
%   * [M_.fname,'.dynamic']
%   * dseries
%   * dsge_likelihood.m
%   * dyn_vech
%   * dynare_resolve
%   * ident_bruteforce
%   * identification_checks
%   * identification_checks_via_subsets
%   * isoctave
%   * get_identification_jacobians (previously getJJ)
%   * matlab_ver_less_than
%   * prior_bounds
%   * set_all_parameters
%   * simulated_moment_uncertainty
%   * stoch_simul
%   * vec
% =========================================================================
% Copyright (C) 2008-2019 Dynare Team
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

global oo_ M_ options_ bayestopt_ estim_params_
persistent ind_J ind_G ind_D ind_dTAU ind_dLRE
% persistence is necessary, because in a MC loop the numerical threshold used may provide vectors of different length, leading to crashes in MC loops

%initialize output structures
ide_hess         = struct(); %asymptotic/simulated information matrix
ide_reducedform  = struct(); %reduced form solution
ide_lre          = struct(); %linear rational expectation model
ide_moments      = struct(); %Iskrev (2010)'s J
ide_spectrum     = struct(); %Qu and Tkachenko (2012)'s G
ide_minimal      = struct(); %Komunjer and Ng (2011)'s D
derivatives_info = struct(); %storage for Jacobians used in dsge_likelihood.m for (analytical or simulated) asymptotic Gradient and Hession of likelihood

totparam_nbr    = length(params);     %number of all parameters to be checked
modparam_nbr    = length(indpmodel);  %number of model parameters to be checked
stderrparam_nbr = length(indpstderr); %number of stderr parameters to be checked
corrparam_nbr   = size(indpcorr,1);   %number of stderr parameters to be checked
indvobs         = bayestopt_.mf2;     % index of observable variables

if ~isempty(estim_params_)
    %estimated_params block is available, so we are able to use set_all_parameters.m
    M_ = set_all_parameters(params,estim_params_,M_);
end

%get options (see dynare_identification.m for description of options)
nlags                  = options_ident.ar;
advanced               = options_ident.advanced;
replic                 = options_ident.replic;
periods                = options_ident.periods;
max_dim_cova_group     = options_ident.max_dim_cova_group;
normalize_jacobians    = options_ident.normalize_jacobians;
tol_deriv              = options_ident.tol_deriv;
tol_rank               = options_ident.tol_rank;
tol_sv                 = options_ident.tol_sv;
no_identification_strength             = options_ident.no_identification_strength;
no_identification_reducedform          = options_ident.no_identification_reducedform;
no_identification_moments              = options_ident.no_identification_moments;
no_identification_minimal              = options_ident.no_identification_minimal;
no_identification_spectrum             = options_ident.no_identification_spectrum;
checks_via_subsets     = options_ident.checks_via_subsets;

[I,~] = find(M_.lead_lag_incidence'); %I is used to select nonzero columns of the Jacobian of endogenous variables in dynamic model files

%Compute linear approximation and matrices A and B of the transition equation for all endogenous variables in DR order
[A, B, ~, info, M_, options_, oo_] = dynare_resolve(M_,options_,oo_);

if info(1) == 0 %no errors in solution
    oo0      = oo_;                         %store oo_ structure
    y0       = oo_.dr.ys(oo_.dr.order_var); %steady state in DR order
    yy0      = oo_.dr.ys(I);                %steady state of dynamic (endogenous and auxiliary variables) in DR order
    ex0      = oo_.exo_steady_state';       %steady state of exogenous variables in declaration order
    param0   = M_.params;                   %parameter values at which to evaluate dynamic, static and param_derivs files
    Sigma_e0 = M_.Sigma_e;                  %store covariance matrix of exogenous shocks
    
    %compute vectorized reduced-form solution in DR order
    tau = [y0; vec(A); dyn_vech(B*Sigma_e0*B')];
    
    %compute vectorized linear rational expectations model in DR order
    [~, g1 ] = feval([M_.fname,'.dynamic'], yy0, repmat(ex0, M_.maximum_exo_lag+M_.maximum_exo_lead+1), param0, oo_.dr.ys, 1);
        %g1 is [endo_nbr by (dynamicvar_nbr+exo_nbr)] first derivative (wrt all endogenous, exogenous and auxiliary variables) of dynamic model equations, i.e. df/d[yy0;ex0], in DR order
    lre = [y0; vec(g1)]; %add steady state in DR order and vectorize

    %Compute Jacobians for identification analysis either analytically or numerically dependent on analytic_derivation_mode in options_ident
    [J, G, D, dTAU, dLRE, dA, dOm, dYss, MOMENTS] = get_identification_jacobians(A, B, estim_params_, M_, oo0, options_, options_ident, indpmodel, indpstderr, indpcorr, indvobs);
    if isempty(D)
        % Komunjer and Ng is not computed if (1) minimality conditions are not fullfilled or (2) there are more shocks and measurement errors than observables, so we need to reset options
        no_identification_minimal = 1;
        options_ident.no_identification_minimal = 1;
    end
    
    %store derivatives for use in dsge_likelihood.m
    derivatives_info.DT   = dA;
    derivatives_info.DOm  = dOm;
    derivatives_info.DYss = dYss;
    
    if init
        %check stationarity
        if ~no_identification_moments
            ind_J = (find(max(abs(J'),[],1) > tol_deriv)); %index for non-zero derivatives
            if isempty(ind_J) && any(any(isnan(J)))
                error('There are NaN in the JJ matrix. Please check whether your model has units roots and you forgot to set diffuse_filter=1.' )
            end
            if any(any(isnan(MOMENTS)))
                error('There are NaN''s in the theoretical moments: make sure that for non-stationary models stationary transformations of non-stationary observables are used for checking identification. [TIP: use first differences].')
            end
        end
        if ~no_identification_spectrum
            ind_G = (find(max(abs(G'),[],1) > tol_deriv)); %index for non-zero derivatives
            if isempty(ind_G) && any(any(isnan(G)))
                warning_QuTkachenko = 'WARNING: There are NaN in Qu and Tkachenko (2012)''s G matrix. Please check whether your model has units roots (diffuse_filter=1 does not work yet for the G matrix).\n';
                warning_QuTkachenko = [warning_QuTkachenko '         Skip identification analysis based on spectrum.\n'];
                fprintf(warning_QuTkachenko);
                %reset options to neither display nor plot Qu and Tkachenko's G anymore
                no_identification_spectrum = 1;
                options_ident.no_identification_spectrum = 1;
            end
        end
        if ~no_identification_minimal
            ind_D = (find(max(abs(D'),[],1) > tol_deriv)); %index for non-zero derivatives
            if isempty(ind_D) && any(any(isnan(D)))
                warning_KomunjerNg = 'WARNING: There are NaN in Komunjer and Ng (2011)''s D matrix. Please check whether your model has units roots and you forgot to set diffuse_filter=1.\n';
                warning_KomunjerNg = [warning_KomunjerNg '         Skip identification analysis based on minimal system.\n'];
                fprintf(warning_KomunjerNg);
                %reset options to neither display nor plot Komunjer and Ng's D anymore
                no_identification_minimal = 1;
                options_ident.no_identification_minimal = 1;
            end
        end
        if no_identification_moments && no_identification_minimal && no_identification_spectrum
            %error if all three criteria fail
            error('identification_analyis: Stationarity condition(s) failed and/or diffuse_filter option missing');
        end
        
        % Check order conditions
        if ~no_identification_moments
            %check order condition of Iskrev (2010)
            while length(ind_J) < totparam_nbr && nlags < 10
                %Try to add lags to autocovariogram if order condition fails
                disp('The number of moments with non-zero derivative is smaller than the number of parameters')
                disp(['Try increasing ar = ', int2str(nlags+1)])
                nlags = nlags + 1;
                options_ident.no_identification_minimal  = 1; %do not recompute D
                options_ident.no_identification_spectrum = 1; %do not recompute G
                options_ident.ar = nlags;     %store new lag number
                options_.ar      = nlags;     %store new lag number
                [J, ~, ~, dTAU, dLRE, dA, dOm, dYss, MOMENTS] = get_identification_jacobians(A, B, estim_params_, M_, oo0, options_, options_ident, indpmodel, indpstderr, indpcorr, indvobs);
                ind_J = (find(max(abs(J'),[],1) > tol_deriv)); %new index with non-zero derivatives
                %store derivatives for use in dsge_likelihood.m
                derivatives_info.DT   = dA;
                derivatives_info.DOm  = dOm;
                derivatives_info.DYss = dYss;                
            end
            options_ident.no_identification_minimal  = no_identification_minimal;  % reset option to original choice
            options_ident.no_identification_spectrum = no_identification_spectrum; % reset option to original choice
            if length(ind_J) < totparam_nbr
                warning_Iskrev = 'WARNING: Order condition for Iskrev (2010) failed: There are not enough moments and too many parameters.\n';
                warning_Iskrev = [warning_Iskrev '         The number of moments with non-zero derivative is smaller than the number of parameters up to 10 lags.\n'];
                warning_Iskrev = [warning_Iskrev '         Either reduce the list of parameters, or further increase ar, or increase number of observables.\n'];
                warning_Iskrev = [warning_Iskrev '         Skip identification analysis based on moments.\n'];
                fprintf(warning_Iskrev);
                %reset options to neither display nor plot Iskrev's J anymore
                no_identification_moments = 1;
                options_ident.no_identification_moments = 1;
            end
        end
        if ~no_identification_minimal
            if length(ind_D) < size(D,2)                
                warning_KomunjerNg = 'WARNING: Order condition for Komunjer and Ng (2011) failed: There are too many parameters or too few observable variables.\n';
                warning_KomunjerNg = [warning_KomunjerNg '         The number of minimal system elements with non-zero derivative is smaller than the number of parameters.\n'];
                warning_KomunjerNg = [warning_KomunjerNg '         Either reduce the list of parameters, or increase number of observables.\n'];
                warning_KomunjerNg = [warning_KomunjerNg '         Skip identification analysis based on minimal state space system.\n'];
                fprintf(warning_KomunjerNg);
                %resest options to neither display nor plot Komunjer and Ng's D anymore
                no_identification_minimal = 1;
                options_ident.no_identification_minimal = 1;
            end
        end
        %There is no order condition for Qu and Tkachenko's G        
        if no_identification_moments && no_identification_minimal && no_identification_spectrum
            %error if all three criteria fail
            error('identification_analyis: Order condition(s) failed');
        end
        if ~no_identification_reducedform
            ind_dTAU = (find(max(abs(dTAU'),[],1) > tol_deriv)); %index with non-zero derivatives
        end
        ind_dLRE = (find(max(abs(dLRE'),[],1) > tol_deriv)); %index with non-zero derivatives
    end
    
    LRE = lre(ind_dLRE);
    si_dLRE = (dLRE(ind_dLRE,:)); %focus only on non-zero derivatives
    if ~no_identification_reducedform
        TAU = tau(ind_dTAU);
        si_dTAU = (dTAU(ind_dTAU,:)); %focus only on non-zero derivatives
    end
    
    if ~no_identification_moments
        MOMENTS = MOMENTS(ind_J);
        si_J   = (J(ind_J,:)); %focus only on non-zero derivatives
        %% compute identification strength based on moments
        if ~no_identification_strength && init %only for initialization of persistent vars
            ide_strength_J        = NaN(1,totparam_nbr); %initialize
            ide_strength_J_prior  = NaN(1,totparam_nbr); %initialize
            ide_uncert_unnormaliz = NaN(1,totparam_nbr); %initialize
            if prior_exist
                offset_prior = 0;
                if ~isempty(estim_params_.var_exo) %stderr parameters come first
                    normaliz_prior_std = bayestopt_.p2(1:estim_params_.nvx)'; % normalize with prior standard deviation
                    offset_prior = offset_prior+estim_params_.nvx+estim_params_.nvn;
                else
                    normaliz_prior_std=[]; %initialize
                end
                if ~isempty(estim_params_.corrx) %corr parameters come second
                    normaliz_prior_std = [normaliz_prior_std bayestopt_.p2(offset_prior+1:offset_prior+estim_params_.ncx)']; % normalize with prior standard deviation
                    offset_prior = offset_prior+estim_params_.ncx+estim_params_.ncn;
                end
                if ~isempty(estim_params_.param_vals) %model parameters come third
                    normaliz_prior_std = [normaliz_prior_std bayestopt_.p2(offset_prior+1:offset_prior+estim_params_.np)']; % normalize with prior standard deviation
                end
            else
                normaliz_prior_std = NaN(1,totparam_nbr); %no prior information available, do not normalize
            end
            
            try
                %try to compute asymptotic Hessian for identification strength analysis based on moments
                % reset some options for faster computations
                options_.irf                     = 0;
                options_.noprint                 = 1;
                options_.order                   = 1;
                options_.SpectralDensity.trigger = 0;
                options_.periods                 = periods+100;
                if options_.kalman_algo > 2
                    options_.kalman_algo         = 1;
                end
                analytic_derivation              = options_.analytic_derivation;
                options_.analytic_derivation     = -2; %this sets asy_Hess=1 in dsge_likelihood.m                
                [info, oo_, options_] = stoch_simul(M_, options_, oo_, options_.varobs);
                dataset_ = dseries(oo_.endo_simul(options_.varobs_id,100+1:end)',dates('1Q1'), options_.varobs); %get information on moments
                derivatives_info.no_DLIK = 1;
                bounds = prior_bounds(bayestopt_, options_.prior_trunc); %reset bounds as lb and ub must only be operational during mode-finding
                [fval, info, cost_flag, DLIK, AHess, ys, trend_coeff, M_, options_, bayestopt_, oo_] = dsge_likelihood(params', dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, bounds, oo_, derivatives_info); %non-used output variables need to be set for octave for some reason
                    %note that for the order of parameters in AHess we have: stderr parameters come first, corr parameters second, model parameters third. the order within these blocks corresponds to the order specified in the estimated_params block
                options_.analytic_derivation = analytic_derivation; %reset option
                AHess = -AHess; %take negative of hessian
                if ischar(tol_rank) %if tolerance is specified as 'robust'
                    tol_rankAhess = max(size(AHess))*eps(norm(AHess));
                else
                    tol_rankAhess = tol_rank;
                end
                if min(eig(AHess))<-tol_rankAhess
                    error('identification_analysis: Analytic Hessian is not positive semi-definite!')
                end
                ide_hess.AHess = AHess; %store asymptotic Hessian
                %normalize asymptotic hessian
                deltaM = sqrt(diag(AHess));
                iflag = any((deltaM.*deltaM)==0); %check if all second-order derivatives wrt to a single parameter are nonzero
                tildaM = AHess./((deltaM)*(deltaM')); %this normalization is for numerical purposes
                if iflag || rank(AHess)>rank(tildaM)
                    [ide_hess.cond, ide_hess.rank, ide_hess.ind0, ide_hess.indno, ide_hess.ino, ide_hess.Mco, ide_hess.Pco] = identification_checks(AHess, 0, tol_rank, tol_sv, totparam_nbr);
                else %use normalized version if possible
                    [ide_hess.cond, ide_hess.rank, ide_hess.ind0, ide_hess.indno, ide_hess.ino, ide_hess.Mco, ide_hess.Pco] = identification_checks(tildaM, 0, tol_rank, tol_sv, totparam_nbr);
                end
                indok = find(max(ide_hess.indno,[],1)==0);
                ide_uncert_unnormaliz(indok) = sqrt(diag(inv(AHess(indok,indok))))';
                ind1 = find(ide_hess.ind0);
                cmm = si_J(:,ind1)*((AHess(ind1,ind1))\si_J(:,ind1)'); %covariance matrix of moments
                temp1 = ((AHess(ind1,ind1))\si_dTAU(:,ind1)');
                diag_chh = sum(si_dTAU(:,ind1)'.*temp1)';
                ind1 = ind1(ind1>stderrparam_nbr+corrparam_nbr);
                clre = si_dLRE(:,ind1-stderrparam_nbr-corrparam_nbr)*((AHess(ind1,ind1))\si_dLRE(:,ind1-stderrparam_nbr-corrparam_nbr)');
                flag_score = 1; %this is only used for a title of a plot in plot_identification.m
            catch
                %Asymptotic Hessian via simulation                
                replic = max([replic, length(ind_J)*3]);
                cmm = simulated_moment_uncertainty(ind_J, periods, replic,options_,M_,oo_); %covariance matrix of moments
                sd = sqrt(diag(cmm));
                cc = cmm./(sd*sd');
                if isoctave || matlab_ver_less_than('8.3')
                    [VV,DD] = eig(cc);
                    %fix for older Matlab versions that do not support computing left eigenvalues, see http://mathworks.com/help/releases/R2012b/matlab/ref/eig.html
                    [WW,~] = eig(cc.');
                    WW = conj(WW);
                else
                    [VV,DD,WW] = eig(cc);
                end
                id = find(diag(DD)>tol_deriv);
                siTMP = si_J./repmat(sd,[1 totparam_nbr]);
                MIM = (siTMP'*VV(:,id))*(DD(id,id)\(WW(:,id)'*siTMP));
                clear siTMP;
                ide_hess.AHess = MIM; %store asymptotic hessian
                %normalize asymptotic hessian
                deltaM = sqrt(diag(MIM));
                iflag = any((deltaM.*deltaM)==0);
                tildaM = MIM./((deltaM)*(deltaM'));
                if iflag || rank(MIM)>rank(tildaM)
                    [ide_hess.cond, ide_hess.rank, ide_hess.ind0, ide_hess.indno, ide_hess.ino, ide_hess.Mco, ide_hess.Pco] = identification_checks(MIM, 0, tol_rank, tol_sv, totparam_nbr);
                else %use normalized version if possible
                    [ide_hess.cond, ide_hess.rank, ide_hess.ind0, ide_hess.indno, ide_hess.ino, ide_hess.Mco, ide_hess.Pco] = identification_checks(tildaM, 0, tol_rank, tol_sv, totparam_nbr);
                end
                indok = find(max(ide_hess.indno,[],1)==0);
                ind1 = find(ide_hess.ind0);
                temp1 = ((MIM(ind1,ind1))\si_dTAU(:,ind1)');
                diag_chh = sum(si_dTAU(:,ind1)'.*temp1)';
                ind1 = ind1(ind1>stderrparam_nbr+corrparam_nbr);
                clre = si_dLRE(:,ind1-stderrparam_nbr-corrparam_nbr)*((MIM(ind1,ind1))\si_dLRE(:,ind1-stderrparam_nbr-corrparam_nbr)');
                if ~isempty(indok)
                    ide_uncert_unnormaliz(indok) = (sqrt(diag(inv(tildaM(indok,indok))))./deltaM(indok))'; %sqrt(diag(inv(MIM(indok,indok))))';
                end
                flag_score = 0; %this is only used for a title of a plot
            end % end of computing sample information matrix for identification strength measure
            
            ide_strength_J(indok) = (1./(ide_uncert_unnormaliz(indok)'./abs(params(indok)')));              %this is s_i in Ratto and Iskrev (2011, p.13)
            ide_strength_J_prior(indok) = (1./(ide_uncert_unnormaliz(indok)'./normaliz_prior_std(indok)')); %this is s_i^{prior} in Ratto and Iskrev (2011, p.14)
            sensitivity_zero_pos = find(isinf(deltaM));
            deltaM_prior = deltaM.*abs(normaliz_prior_std'); %this is \Delta_i^{prior} in Ratto and Iskrev (2011, p.14)
            deltaM = deltaM.*abs(params');                   %this is \Delta_i in Ratto and Iskrev (2011, p.14)
            quant = si_J./repmat(sqrt(diag(cmm)),1,totparam_nbr);
            if size(quant,1)==1
                si_Jnorm = abs(quant).*normaliz_prior_std;
            else
                si_Jnorm = vnorm(quant).*normaliz_prior_std;
            end
            iy = find(diag_chh);
            ind_dTAU = ind_dTAU(iy);
            si_dTAU = si_dTAU(iy,:);
            if ~isempty(iy)
                quant = si_dTAU./repmat(sqrt(diag_chh(iy)),1,totparam_nbr);
                if size(quant,1)==1
                    si_dTAUnorm = abs(quant).*normaliz_prior_std;
                else
                    si_dTAUnorm = vnorm(quant).*normaliz_prior_std;
                end
            else
                si_dTAUnorm = [];
            end
            diag_clre = diag(clre);
            iy = find(diag_clre);
            ind_dLRE = ind_dLRE(iy);
            si_dLRE = si_dLRE(iy,:);
            if ~isempty(iy)
                quant = si_dLRE./repmat(sqrt(diag_clre(iy)),1,modparam_nbr);
                if size(quant,1)==1
                    si_dLREnorm = abs(quant).*normaliz_prior_std(stderrparam_nbr+corrparam_nbr+1:end);
                else
                    si_dLREnorm = vnorm(quant).*normaliz_prior_std(stderrparam_nbr+corrparam_nbr+1:end);
                end
            else
                si_dLREnorm=[];
            end
            %store results of identification strength
            ide_hess.ide_strength_J               = ide_strength_J;
            ide_hess.ide_strength_J_prior         = ide_strength_J_prior;
            ide_hess.deltaM                       = deltaM;
            ide_hess.deltaM_prior                 = deltaM_prior;
            ide_hess.sensitivity_zero_pos         = sensitivity_zero_pos;
            ide_hess.identified_parameter_indices = indok;
            ide_hess.flag_score                   = flag_score;            
            
            ide_lre.si_dLREnorm                   = si_dLREnorm;
            ide_moments.si_Jnorm                  = si_Jnorm;            
            ide_reducedform.si_dTAUnorm           = si_dTAUnorm;            
        end %end of identification strength analysis
    end
    
    %% Theoretical identification analysis based on linear rational expectations model, i.e. steady state and Jacobian of dynamic model equations
    if normalize_jacobians
        norm_dLRE = max(abs(si_dLRE),[],2);
        norm_dLRE = norm_dLRE(:,ones(size(dLRE,2),1));
    else
        norm_dLRE = 1;
    end
    % store into structure (not everything is really needed)
    ide_lre.ind_dLRE  = ind_dLRE;
    ide_lre.norm_dLRE = norm_dLRE;
    ide_lre.si_dLRE   = si_dLRE;
    ide_lre.dLRE      = dLRE;
    ide_lre.LRE       = LRE;
    
    %% Theoretical identification analysis based on steady state and reduced-form-solution
    if ~no_identification_reducedform
        if normalize_jacobians
            norm_dTAU = max(abs(si_dTAU),[],2);
            norm_dTAU = norm_dTAU(:,ones(totparam_nbr,1));
        else
            norm_dTAU = 1;
        end
        % store into structure (not everything is really needed)
        ide_reducedform.ind_dTAU  = ind_dTAU;
        ide_reducedform.norm_dTAU = norm_dTAU;
        ide_reducedform.si_dTAU   = si_dTAU;
        ide_reducedform.dTAU      = dTAU;
        ide_reducedform.TAU       = TAU;
    end    
    
    %% Theoretical identification analysis based on Iskrev (2010)'s method, i.e. first two moments
    if ~no_identification_moments
        if normalize_jacobians
            norm_J = max(abs(si_J),[],2);
            norm_J = norm_J(:,ones(totparam_nbr,1));
        else
            norm_J = 1;
        end
        % store into structure (not everything is really needed)
        ide_moments.ind_J     = ind_J;
        ide_moments.norm_J    = norm_J;
        ide_moments.si_J      = si_J;
        ide_moments.J         = J;
        ide_moments.MOMENTS   = MOMENTS;
        
        if advanced
            % here we do not normalize (i.e. set normJ=1) as the OLS in ident_bruteforce is very sensitive to normJ
            [ide_moments.pars, ide_moments.cosnJ] = ident_bruteforce(J(ind_J,:), max_dim_cova_group, options_.TeX, options_ident.name_tex, options_ident.tittxt, options_ident.tol_deriv);
        end
        %here idea is to have the unnormalized S and V, which is then used in plot_identification.m and for prior_mc > 1
        [~, S, V] = svd(J(ind_J,:),0);
        S = diag(S);
        S = [S;zeros(size(J,2)-length(ind_J),1)];
        if totparam_nbr > 8
            ide_moments.S = S([1:4, end-3:end]);
            ide_moments.V = V(:,[1:4, end-3:end]);
        else
            ide_moments.S = S;
            ide_moments.V = V;
        end
    end
    
    %% Theoretical identification analysis based on Komunjer and Ng (2011)'s method, i.e. steady state and observational equivalent spectral densities within a minimal system
    if ~no_identification_minimal
        if normalize_jacobians
            norm_D = max(abs(D(ind_D,:)),[],2);
            norm_D = norm_D(:,ones(size(D,2),1));
        else
            norm_D = 1;
        end
        % store into structure
        ide_minimal.ind_D  = ind_D;
        ide_minimal.norm_D = norm_D;
        ide_minimal.D      = D;
    end
    
    %% Theoretical identification analysis based on Qu and Tkachenko (2012)'s method, i.e. steady state and spectral density
    if ~no_identification_spectrum
        if normalize_jacobians
            tilda_G = zeros(size(G));
            delta_G = sqrt(diag(G(ind_G,ind_G)));
            tilda_G(ind_G,ind_G) = G(ind_G,ind_G)./((delta_G)*(delta_G'));
            norm_G = max(abs(G(ind_G,:)),[],2);
            norm_G = norm_G(:,ones(size(G,2),1));
        else
            tilda_G = G;
            norm_G = 1;
        end
        % store into structure
        ide_spectrum.ind_G  = ind_G;
        ide_spectrum.tilda_G = tilda_G;
        ide_spectrum.G      = G;
        ide_spectrum.norm_G = norm_G;
    end
    
    %% Perform identification checks, i.e. find out which parameters are involved
    if checks_via_subsets
        % identification_checks_via_subsets is only for debugging
        [ide_lre, ide_reducedform, ide_moments, ide_spectrum, ide_minimal] = ...
            identification_checks_via_subsets(ide_lre, ide_reducedform, ide_moments, ide_spectrum, ide_minimal, totparam_nbr, modparam_nbr, options_ident);
    else
        [ide_lre.cond, ide_lre.rank, ide_lre.ind0, ide_lre.indno, ide_lre.ino, ide_lre.Mco, ide_lre.Pco, ide_lre.jweak, ide_lre.jweak_pair] = ...
            identification_checks(dLRE(ind_dLRE,:)./norm_dLRE, 1, tol_rank, tol_sv, modparam_nbr);
        if ~no_identification_reducedform
            [ide_reducedform.cond, ide_reducedform.rank, ide_reducedform.ind0, ide_reducedform.indno, ide_reducedform.ino, ide_reducedform.Mco, ide_reducedform.Pco, ide_reducedform.jweak, ide_reducedform.jweak_pair] = ...
                identification_checks(dTAU(ind_dTAU,:)./norm_dTAU, 1, tol_rank, tol_sv, totparam_nbr);
        end
        if ~no_identification_moments
            [ide_moments.cond, ide_moments.rank, ide_moments.ind0, ide_moments.indno, ide_moments.ino, ide_moments.Mco, ide_moments.Pco, ide_moments.jweak, ide_moments.jweak_pair] = ...
                identification_checks(J(ind_J,:)./norm_J, 1, tol_rank, tol_sv, totparam_nbr);
        end                
        if ~no_identification_minimal
            [ide_minimal.cond, ide_minimal.rank, ide_minimal.ind0, ide_minimal.indno, ide_minimal.ino, ide_minimal.Mco, ide_minimal.Pco, ide_minimal.jweak, ide_minimal.jweak_pair] = ...
                identification_checks(D(ind_D,:)./norm_D, 2, tol_rank, tol_sv, totparam_nbr);
        end        
        if ~no_identification_spectrum            
            [ide_spectrum.cond, ide_spectrum.rank, ide_spectrum.ind0, ide_spectrum.indno, ide_spectrum.ino, ide_spectrum.Mco, ide_spectrum.Pco, ide_spectrum.jweak, ide_spectrum.jweak_pair] = ...
                identification_checks(tilda_G, 3, tol_rank, tol_sv, totparam_nbr);            
        end
    end
end