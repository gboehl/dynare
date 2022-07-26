function [ide_moments, ide_spectrum, ide_minimal, ide_hess, ide_reducedform, ide_dynamic, derivatives_info, info, error_indicator] = identification_analysis(params, indpmodel, indpstderr, indpcorr, options_ident, dataset_info, prior_exist, init)
% [ide_moments, ide_spectrum, ide_minimal, ide_hess, ide_reducedform, ide_dynamic, derivatives_info, info, error_indicator] = identification_analysis(params, indpmodel, indpstderr, indpcorr, options_ident, dataset_info, prior_exist, init)
% -------------------------------------------------------------------------
% This function wraps all identification analysis, i.e. it
% (1) wraps functions for the theoretical identification analysis based on moments (Iskrev, 2010),
% spectrum (Qu and Tkachenko, 2012), minimal system (Komunjer and Ng, 2011),  information matrix,
% reduced-form solution and dynamic model derivatives (Ratto and Iskrev, 2011).
% (2) computes the identification strength based on moments (Ratto and Iskrev, 2011)
% (3) checks which parameters are involved.
% If options_ident.order>1, then the identification analysis is based on
% Mutschler (2015), i.e. the pruned state space system and the corresponding
% moments, spectrum, reduced-form solution and dynamic model derivatives
% =========================================================================
% INPUTS
%    * params             [mc_sample_nbr by totparam_nbr]
%                         parameter values for identification checks
%    * indpmodel          [modparam_nbr by 1]
%                         index of model parameters for which identification is checked
%    * indpstderr         [stderrparam_nbr by 1]
%                         index of stderr parameters for which identification is checked
%    * indpcorr           [corrparam_nbr by 2]
%                         matrix of corr parmeters for which identification is checked
%    * options_ident      [structure]
%                         identification options
%    * dataset_info       [structure]
%                         various information about the dataset (descriptive statistics and missing observations) for Kalman Filter
%    * prior_exist        [integer]
%                         1: prior exists. Identification is checked for stderr, corr and model parameters as declared in estimated_params block
%                         0: prior does not exist. Identification is checked for all stderr and model parameters, but no corr parameters
%    * init               [integer]
%                         flag for initialization of persistent vars. This is needed if one may want to make more calls to identification in the same mod file
%    * error_indicator    [structure] 
%                         boolean information on errors (1 is an error, 0 is no error) while computing the criteria stored in fields identification_strength, identification_reducedform, identification_moments, identification_minimal, identification_spectrum
% -------------------------------------------------------------------------
% OUTPUTS
%    * ide_moments        [structure]
%                         identification results using theoretical moments (Iskrev, 2010; Mutschler, 2015)
%    * ide_spectrum       [structure]
%                         identification results using spectrum (Qu and Tkachenko, 2012; Mutschler, 2015)
%    * ide_minimal        [structure]
%                         identification results using theoretical mean and minimal system (Komunjer and Ng, 2011)
%    * ide_hess           [structure]
%                         identification results using asymptotic Hessian (Ratto and Iskrev, 2011)
%    * ide_reducedform    [structure]
%                         identification results using steady state and reduced form solution (Ratto and Iskrev, 2011)
%    * ide_dynamic        [structure]
%                         identification results using steady state and dynamic model derivatives (Ratto and Iskrev, 2011)
%    * derivatives_info   [structure]
%                         info about first-order perturbation derivatives, used in dsge_likelihood.m
%    * info               [integer]
%                         output from dynare_resolve
%    * error_indicator    [structure]
%                         indicator on problems
% -------------------------------------------------------------------------
% This function is called by
%   * dynare_identification.m
% -------------------------------------------------------------------------
% This function calls
%   * [M_.fname,'.dynamic']
%   * dseries
%   * dsge_likelihood.m
%   * dyn_vech
%   * ident_bruteforce
%   * identification_checks
%   * identification_checks_via_subsets
%   * isoctave
%   * get_identification_jacobians (previously getJJ)
%   * matlab_ver_less_than
%   * prior_bounds
%   * resol
%   * set_all_parameters
%   * simulated_moment_uncertainty
%   * stoch_simul
%   * vec
% =========================================================================
% Copyright (C) 2008-2021 Dynare Team
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

global oo_ M_ options_ bayestopt_ estim_params_
persistent ind_dMOMENTS ind_dREDUCEDFORM ind_dDYNAMIC
% persistent indices are necessary, because in a MC loop the numerical threshold
% used may provide vectors of different length, leading to crashes in MC loops

%initialize output structures
ide_hess         = struct(); %Identification structure based on asymptotic/simulated information matrix
ide_reducedform  = struct(); %Identification structure based on steady state and reduced form solution
ide_dynamic      = struct(); %Identification structure based on steady state and dynamic model derivatives
ide_moments      = struct(); %Identification structure based on first two moments (Iskrev, 2010; Mutschler, 2015)
ide_spectrum     = struct(); %Identification structure based on Gram matrix of Jacobian of spectral density plus Gram matrix of Jacobian of steady state (Qu and Tkachenko, 2012; Mutschler, 2015)
ide_minimal      = struct(); %Identification structure based on mean and minimal system (Komunjer and Ng, 2011)
derivatives_info = struct(); %storage for first-order perturbation Jacobians used in dsge_likelihood.m

totparam_nbr    = length(params);     %number of all parameters to be checked
modparam_nbr    = length(indpmodel);  %number of model parameters to be checked
stderrparam_nbr = length(indpstderr); %number of stderr parameters to be checked
corrparam_nbr   = size(indpcorr,1);   %number of stderr parameters to be checked
indvobs         = bayestopt_.mf2;     %index of observable variables

if ~isempty(estim_params_)
    %estimated_params block is available, so we are able to use set_all_parameters.m
    M_ = set_all_parameters(params,estim_params_,M_);
end

%get options (see dynare_identification.m for description of options)
nlags               = options_ident.ar;
advanced            = options_ident.advanced;
replic              = options_ident.replic;
periods             = options_ident.periods;
max_dim_cova_group  = options_ident.max_dim_cova_group;
normalize_jacobians = options_ident.normalize_jacobians;
checks_via_subsets  = options_ident.checks_via_subsets;
tol_deriv           = options_ident.tol_deriv;
tol_rank            = options_ident.tol_rank;
tol_sv              = options_ident.tol_sv;
error_indicator.identification_strength=0;
error_indicator.identification_reducedform=0;
error_indicator.identification_moments=0;
error_indicator.identification_minimal=0;
error_indicator.identification_spectrum=0;

%Compute linear approximation and fill dr structure
[oo_.dr,info,M_,oo_] = compute_decision_rules(M_,options_,oo_);

if info(1) == 0 %no errors in solution
    % Compute parameter Jacobians for identification analysis
    [MEAN, dMEAN, REDUCEDFORM, dREDUCEDFORM, DYNAMIC, dDYNAMIC, MOMENTS, dMOMENTS, dSPECTRUM, dSPECTRUM_NO_MEAN, dMINIMAL, derivatives_info] = get_identification_jacobians(estim_params_, M_, oo_, options_, options_ident, indpmodel, indpstderr, indpcorr, indvobs);
    if isempty(dMINIMAL)
        % Komunjer and Ng is not computed if (1) minimality conditions are not fullfilled or (2) there are more shocks and measurement errors than observables, so we need to reset options
        error_indicator.identification_minimal = 1;
        %options_ident.no_identification_minimal = 1;
    end

    if init
        %check stationarity
        if ~options_ident.no_identification_moments
            if any(any(isnan(MOMENTS)))
                if options_.diffuse_filter == 1 % use options_ as it inherits diffuse_filter from options_ident if set by user
                    error('There are NaN''s in the theoretical moments. Make sure that for non-stationary models stationary transformations of non-stationary observables are used for checking identification. [TIP: use first differences].')
                else
                    error('There are NaN''s in the theoretical moments. Please check whether your model has units roots, and you forgot to set diffuse_filter=1.' )
                end
                error_indicator.identification_moments=1;
            end
            ind_dMOMENTS = (find(max(abs(dMOMENTS'),[],1) > tol_deriv)); %index for non-zero rows
            if isempty(ind_dMOMENTS) && any(any(isnan(dMOMENTS)))                
                error('There are NaN in the dMOMENTS matrix.')
                error_indicator.identification_moments=1;
            end            
        end
        if ~options_ident.no_identification_spectrum
            ind_dSPECTRUM = (find(max(abs(dSPECTRUM'),[],1) > tol_deriv)); %index for non-zero rows
            if any(any(isnan(dSPECTRUM)))
                warning_SPECTRUM = 'WARNING: There are NaN in the dSPECTRUM matrix. Note that identification based on spectrum does not support non-stationary models (yet).\n';
                warning_SPECTRUM = [warning_SPECTRUM '         Skip identification analysis based on spectrum.\n'];
                fprintf(warning_SPECTRUM);
                %set indicator to neither display nor plot dSPECTRUM anymore
                error_indicator.identification_spectrum = 1;                
            end
        end
        if ~options_ident.no_identification_minimal
            ind_dMINIMAL = (find(max(abs(dMINIMAL'),[],1) > tol_deriv)); %index for non-zero rows
            if any(any(isnan(dMINIMAL)))
                warning_MINIMAL = 'WARNING: There are NaN in the dMINIMAL matrix. Note that identification based on minimal system does not support non-stationary models (yet).\n';
                warning_MINIMAL = [warning_MINIMAL '         Skip identification analysis based on minimal system.\n'];
                fprintf(warning_MINIMAL);
                %set indicator to neither display nor plot dMINIMAL anymore
                error_indicator.identification_minimal = 1;
            end
        end
        %The following cannot be reached yet due to erroring out when
        %error_indicator.identification_moments is triggered
        if error_indicator.identification_moments && error_indicator.identification_minimal && error_indicator.identification_spectrum
            %display error if all three criteria fail            
            error(sprintf('identification_analyis: Stationarity condition(s) failed and/or diffuse_filter option missing.\nMake sure that for non-stationary models stationary transformations of non-stationary observables are used for checking identification.\n[TIP: use first differences].'));
        end

        % Check order conditions
        if ~options_ident.no_identification_moments && ~error_indicator.identification_moments
            %check order condition of Iskrev (2010)
            while length(ind_dMOMENTS) < totparam_nbr && nlags < 10
                %Try to add lags to autocovariogram if order condition fails
                disp('The number of moments with non-zero derivative is smaller than the number of parameters')
                disp(['Try increasing ar = ', int2str(nlags+1)])
                nlags = nlags + 1;
                options_ident_local=options_ident;
                options_ident_local.no_identification_minimal  = 1; %do not recompute dMINIMAL
                options_ident_local.no_identification_spectrum = 1; %do not recompute dSPECTRUM
                options_ident_local.ar = nlags;     %store new lag number
                options_.ar      = nlags;           %store new lag number
                [~, ~, ~, ~, ~, ~, MOMENTS, dMOMENTS, ~, ~, ~, ~] = get_identification_jacobians(estim_params_, M_, oo_, options_, options_ident_local, indpmodel, indpstderr, indpcorr, indvobs);

                ind_dMOMENTS = (find(max(abs(dMOMENTS'),[],1) > tol_deriv)); %new index with non-zero rows
            end
            if length(ind_dMOMENTS) < totparam_nbr
                warning_MOMENTS = 'WARNING: Order condition for dMOMENTS failed: There are not enough moments and too many parameters.\n';
                warning_MOMENTS = [warning_MOMENTS '         The number of moments with non-zero derivative is smaller than the number of parameters up to 10 lags.\n'];
                warning_MOMENTS = [warning_MOMENTS '         Either reduce the list of parameters, or further increase ar, or increase number of observables.\n'];
                warning_MOMENTS = [warning_MOMENTS '         Skip identification analysis based on moments.\n'];
                warning_MOMENTS = [warning_MOMENTS '         Skip identification strength analysis.\n'];
                fprintf(warning_MOMENTS);
                %set indicator to neither display nor plot dMOMENTS anymore
                error_indicator.identification_moments = 1;
                %options_ident.no_identification_moments = 1;
                error_indicator.identification_strength = 1;
                %options_ident.no_identification_strength = 1;
            end
        end
        if ~options_ident.no_identification_minimal && ~error_indicator.identification_minimal
            if length(ind_dMINIMAL) < size(dMINIMAL,2)
                warning_MINIMAL = 'WARNING: Order condition for dMINIMAL failed: There are too many parameters or too few observable variables.\n';
                warning_MINIMAL = [warning_MINIMAL '         The number of minimal system elements with non-zero derivative is smaller than the number of parameters.\n'];
                warning_MINIMAL = [warning_MINIMAL '         Either reduce the list of parameters, or increase number of observables.\n'];
                warning_MINIMAL = [warning_MINIMAL '         Skip identification analysis based on minimal state space system.\n'];
                fprintf(warning_MINIMAL);
                %set indicator to neither display nor plot dMINIMAL anymore
                error_indicator.identification_minimal = 1;                
            end
        end
        %Note that there is no order condition for dSPECTRUM, as the matrix is always of dimension totparam_nbr by totparam_nbr
        if error_indicator.identification_moments && error_indicator.identification_minimal && error_indicator.identification_spectrum
            %error if all three criteria fail
            error('identification_analyis: Order condition(s) failed');
        end
        if ~options_ident.no_identification_reducedform && ~error_indicator.identification_reducedform
            ind_dREDUCEDFORM = (find(max(abs(dREDUCEDFORM'),[],1) > tol_deriv)); %index with non-zero rows
        end
        ind_dDYNAMIC = (find(max(abs(dDYNAMIC'),[],1) > tol_deriv)); %index with non-zero rows
    end

    DYNAMIC = DYNAMIC(ind_dDYNAMIC); %focus only on non-zero entries
    si_dDYNAMIC = (dDYNAMIC(ind_dDYNAMIC,:)); %focus only on non-zero rows
    if ~options_ident.no_identification_reducedform && ~error_indicator.identification_reducedform
        REDUCEDFORM = REDUCEDFORM(ind_dREDUCEDFORM); %focus only on non-zero entries
        si_dREDUCEDFORM = (dREDUCEDFORM(ind_dREDUCEDFORM,:)); %focus only on non-zero rows
    end

    if ~options_ident.no_identification_moments && ~error_indicator.identification_moments
        MOMENTS = MOMENTS(ind_dMOMENTS); %focus only on non-zero entries
        si_dMOMENTS   = (dMOMENTS(ind_dMOMENTS,:)); %focus only on non-zero derivatives
        %% MOMENTS IDENTIFICATION STRENGTH ANALYSIS
        if ~options_ident.no_identification_strength && ~error_indicator.identification_strength && init %only for initialization of persistent vars
            ide_strength_dMOMENTS        = NaN(1,totparam_nbr); %initialize
            ide_strength_dMOMENTS_prior  = NaN(1,totparam_nbr); %initialize
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
                if options_.order > 1
                    error('IDENTIFICATION STRENGTH: Analytic computation of Hessian is not available for ''order>1''. Identification strength is based on simulated moment uncertainty');
                end                
                % reset some options for faster computations
                options_.irf                     = 0;
                options_.noprint                 = 1;
                options_.SpectralDensity.trigger = 0;
                options_.periods                 = periods+100;
                if options_.kalman_algo > 2
                    options_.kalman_algo         = 1;
                end
                analytic_derivation              = options_.analytic_derivation;
                options_.analytic_derivation     = -2; %this sets asy_Hess=1 in dsge_likelihood.m                
                [info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, options_.varobs);
                dataset_ = dseries(oo_.endo_simul(options_.varobs_id,100+1:end)',dates('1Q1'), options_.varobs); %get information on moments
                derivatives_info.no_DLIK = 1;
                bounds = prior_bounds(bayestopt_, options_.prior_trunc); %reset bounds as lb and ub must only be operational during mode-finding 
                %note that for order>1 we do not provide any information on DT,DYss,DOM in derivatives_info, such that dsge_likelihood creates an error. Therefore the computation will be based on simulated_moment_uncertainty for order>1.
                [fval, info, cost_flag, DLIK, AHess, ys, trend_coeff, M_, options_, bayestopt_, oo_] = dsge_likelihood(params', dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, bounds, oo_, derivatives_info); %non-used output variables need to be set for octave for some reason
                    %note that for the order of parameters in AHess we have: stderr parameters come first, corr parameters second, model parameters third. the order within these blocks corresponds to the order specified in the estimated_params block
                options_.analytic_derivation = analytic_derivation; %reset option
                AHess = -AHess; %take negative of hessian
                if min(eig(AHess))<-tol_rank
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
                cmm = si_dMOMENTS(:,ind1)*((AHess(ind1,ind1))\si_dMOMENTS(:,ind1)'); %covariance matrix of moments
                temp1 = ((AHess(ind1,ind1))\si_dREDUCEDFORM(:,ind1)');
                diag_chh = sum(si_dREDUCEDFORM(:,ind1)'.*temp1)';
                ind1 = ind1(ind1>stderrparam_nbr+corrparam_nbr);
                cdynamic = si_dDYNAMIC(:,ind1-stderrparam_nbr-corrparam_nbr)*((AHess(ind1,ind1))\si_dDYNAMIC(:,ind1-stderrparam_nbr-corrparam_nbr)');
                flag_score = 1; %this is used for the title in plot_identification.m
            catch
                %Asymptotic Hessian via simulation
                if options_.order > 1       
                    % reset some options for faster computations
                    options_.irf                     = 0;
                    options_.noprint                 = 1;
                    options_.SpectralDensity.trigger = 0;
                    options_.periods                 = periods+100;
                end                
                replic = max([replic, length(ind_dMOMENTS)*3]);
                cmm = simulated_moment_uncertainty(ind_dMOMENTS, periods, replic,options_,M_,oo_); %covariance matrix of moments
                sd = sqrt(diag(cmm));
                cc = cmm./(sd*sd');
                if isoctave
                    [VV,DD] = eig(cc);
                    %fix for older Matlab versions that do not support computing left eigenvalues, see http://mathworks.com/help/releases/R2012b/matlab/ref/eig.html
                    [WW,~] = eig(cc.');
                    WW = conj(WW);
                else
                    [VV,DD,WW] = eig(cc);
                end
                id = find(diag(DD)>tol_deriv);
                siTMP = si_dMOMENTS./repmat(sd,[1 totparam_nbr]);
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
                temp1 = ((MIM(ind1,ind1))\si_dREDUCEDFORM(:,ind1)');
                diag_chh = sum(si_dREDUCEDFORM(:,ind1)'.*temp1)';
                ind1 = ind1(ind1>stderrparam_nbr+corrparam_nbr);
                cdynamic = si_dDYNAMIC(:,ind1-stderrparam_nbr-corrparam_nbr)*((MIM(ind1,ind1))\si_dDYNAMIC(:,ind1-stderrparam_nbr-corrparam_nbr)');
                if ~isempty(indok)
                    ide_uncert_unnormaliz(indok) = (sqrt(diag(inv(tildaM(indok,indok))))./deltaM(indok))'; %sqrt(diag(inv(MIM(indok,indok))))';
                end
                flag_score = 0; %this is used for the title in plot_identification.m
            end % end of computing sample information matrix for identification strength measure

            ide_strength_dMOMENTS(indok) = (1./(ide_uncert_unnormaliz(indok)'./abs(params(indok)')));              %this is s_i in Ratto and Iskrev (2011, p.13)
            ide_strength_dMOMENTS_prior(indok) = (1./(ide_uncert_unnormaliz(indok)'./normaliz_prior_std(indok)')); %this is s_i^{prior} in Ratto and Iskrev (2011, p.14)
            sensitivity_zero_pos = find(isinf(deltaM));
            deltaM_prior = deltaM.*abs(normaliz_prior_std'); %this is \Delta_i^{prior} in Ratto and Iskrev (2011, p.14)
            deltaM = deltaM.*abs(params');                   %this is \Delta_i in Ratto and Iskrev (2011, p.14)
            quant = si_dMOMENTS./repmat(sqrt(diag(cmm)),1,totparam_nbr);
            if size(quant,1)==1
                si_dMOMENTSnorm = abs(quant).*normaliz_prior_std;
            else
                si_dMOMENTSnorm = vnorm(quant).*normaliz_prior_std;
            end
            iy = find(diag_chh);
            ind_dREDUCEDFORM = ind_dREDUCEDFORM(iy);
            si_dREDUCEDFORM = si_dREDUCEDFORM(iy,:);
            if ~isempty(iy)
                quant = si_dREDUCEDFORM./repmat(sqrt(diag_chh(iy)),1,totparam_nbr);
                if size(quant,1)==1
                    si_dREDUCEDFORMnorm = abs(quant).*normaliz_prior_std;
                else
                    si_dREDUCEDFORMnorm = vnorm(quant).*normaliz_prior_std;
                end
            else
                si_dREDUCEDFORMnorm = [];
            end
            diag_cdynamic = diag(cdynamic);
            iy = find(diag_cdynamic);
            ind_dDYNAMIC = ind_dDYNAMIC(iy);
            si_dDYNAMIC = si_dDYNAMIC(iy,:);
            if ~isempty(iy)
                quant = si_dDYNAMIC./repmat(sqrt(diag_cdynamic(iy)),1,modparam_nbr);
                if size(quant,1)==1
                    si_dDYNAMICnorm = abs(quant).*normaliz_prior_std(stderrparam_nbr+corrparam_nbr+1:end);
                else
                    si_dDYNAMICnorm = vnorm(quant).*normaliz_prior_std(stderrparam_nbr+corrparam_nbr+1:end);
                end
            else
                si_dDYNAMICnorm=[];
            end
            %store results of identification strength
            ide_hess.ide_strength_dMOMENTS        = ide_strength_dMOMENTS;
            ide_hess.ide_strength_dMOMENTS_prior  = ide_strength_dMOMENTS_prior;
            ide_hess.deltaM                       = deltaM;
            ide_hess.deltaM_prior                 = deltaM_prior;
            ide_hess.sensitivity_zero_pos         = sensitivity_zero_pos;
            ide_hess.identified_parameter_indices = indok;
            ide_hess.flag_score                   = flag_score;            
            
            ide_dynamic.si_dDYNAMICnorm           = si_dDYNAMICnorm;
            ide_moments.si_dMOMENTSnorm           = si_dMOMENTSnorm;            
            ide_reducedform.si_dREDUCEDFORMnorm   = si_dREDUCEDFORMnorm;            
        end %end of identification strength analysis
    end
    
%% Normalization of Jacobians 
% For Dynamic, ReducedForm, Moment and Minimal Jacobian: rescale each row by its largest element in absolute value
% For Spectrum: transform into correlation-type matrices (as above with AHess)
    if normalize_jacobians
        norm_dDYNAMIC = max(abs(si_dDYNAMIC),[],2);
        norm_dDYNAMIC = norm_dDYNAMIC(:,ones(size(dDYNAMIC,2),1));
    else
        norm_dDYNAMIC = 1;
    end
    % store into structure (not everything is used later on)
    ide_dynamic.ind_dDYNAMIC  = ind_dDYNAMIC;
    ide_dynamic.norm_dDYNAMIC = norm_dDYNAMIC;
    ide_dynamic.si_dDYNAMIC   = si_dDYNAMIC;
    ide_dynamic.dDYNAMIC      = dDYNAMIC;
    ide_dynamic.DYNAMIC       = DYNAMIC;

    if ~options_ident.no_identification_reducedform && ~error_indicator.identification_reducedform
        if normalize_jacobians
            norm_dREDUCEDFORM = max(abs(si_dREDUCEDFORM),[],2);
            norm_dREDUCEDFORM = norm_dREDUCEDFORM(:,ones(totparam_nbr,1));
        else
            norm_dREDUCEDFORM = 1;
        end
        % store into structure (not everything is used later on)
        ide_reducedform.ind_dREDUCEDFORM  = ind_dREDUCEDFORM;
        ide_reducedform.norm_dREDUCEDFORM = norm_dREDUCEDFORM;
        ide_reducedform.si_dREDUCEDFORM   = si_dREDUCEDFORM;
        ide_reducedform.dREDUCEDFORM      = dREDUCEDFORM;
        ide_reducedform.REDUCEDFORM       = REDUCEDFORM;
    end

    if ~options_ident.no_identification_moments && ~error_indicator.identification_moments
        if normalize_jacobians
            norm_dMOMENTS = max(abs(si_dMOMENTS),[],2);
            norm_dMOMENTS = norm_dMOMENTS(:,ones(totparam_nbr,1));
        else
            norm_dMOMENTS = 1;
        end
        % store into structure (not everything is used later on)
        ide_moments.ind_dMOMENTS  = ind_dMOMENTS;
        ide_moments.norm_dMOMENTS = norm_dMOMENTS;
        ide_moments.si_dMOMENTS   = si_dMOMENTS;
        ide_moments.dMOMENTS      = dMOMENTS;
        ide_moments.MOMENTS       = MOMENTS;

        if advanced
            % here we do not normalize (i.e. we set norm_dMOMENTS=1) as the OLS in ident_bruteforce is very sensitive to norm_dMOMENTS
            [ide_moments.pars, ide_moments.cosndMOMENTS] = ident_bruteforce(dMOMENTS(ind_dMOMENTS,:), max_dim_cova_group, options_.TeX, options_ident.name_tex, options_ident.tittxt, tol_deriv);
        end

        %here we focus on the unnormalized S and V, which is then used in plot_identification.m and for prior_mc > 1
        [~, S, V] = svd(dMOMENTS(ind_dMOMENTS,:),0);
        if size(S,1) == 1
            S = S(1); % edge case that S is not a matrix but a row vector
        else
            S = diag(S);
        end
        S = [S;zeros(size(dMOMENTS,2)-length(ind_dMOMENTS),1)];
        if totparam_nbr > 8
            ide_moments.S = S([1:4, end-3:end]);
            ide_moments.V = V(:,[1:4, end-3:end]);
        else
            ide_moments.S = S;
            ide_moments.V = V;
        end
    end

    if ~options_ident.no_identification_minimal && ~error_indicator.identification_minimal
        if normalize_jacobians
            ind_dMINIMAL = (find(max(abs(dMINIMAL'),[],1) > tol_deriv)); %index for non-zero rows
            norm_dMINIMAL = max(abs(dMINIMAL(ind_dMINIMAL,:)),[],2);
            norm_dMINIMAL = norm_dMINIMAL(:,ones(size(dMINIMAL,2),1));
        else
            norm_dMINIMAL = 1;
        end
        % store into structure (not everything is used later on)
        ide_minimal.ind_dMINIMAL  = ind_dMINIMAL;
        ide_minimal.norm_dMINIMAL = norm_dMINIMAL;
        ide_minimal.dMINIMAL      = dMINIMAL;
    end

    if ~options_ident.no_identification_spectrum && ~error_indicator.identification_spectrum
        if normalize_jacobians
            ind_dSPECTRUM = (find(max(abs(dSPECTRUM'),[],1) > tol_deriv)); %index for non-zero rows
            tilda_dSPECTRUM = zeros(size(dSPECTRUM));
            delta_dSPECTRUM = sqrt(diag(dSPECTRUM(ind_dSPECTRUM,ind_dSPECTRUM)));
            tilda_dSPECTRUM(ind_dSPECTRUM,ind_dSPECTRUM) = dSPECTRUM(ind_dSPECTRUM,ind_dSPECTRUM)./((delta_dSPECTRUM)*(delta_dSPECTRUM'));
            norm_dSPECTRUM = max(abs(dSPECTRUM(ind_dSPECTRUM,:)),[],2);
            norm_dSPECTRUM = norm_dSPECTRUM(:,ones(size(dSPECTRUM,2),1));
        else
            tilda_dSPECTRUM = dSPECTRUM;
            norm_dSPECTRUM = 1;
        end
        % store into structure (not everything is used later on)
        ide_spectrum.ind_dSPECTRUM   = ind_dSPECTRUM;
        ide_spectrum.norm_dSPECTRUM  = norm_dSPECTRUM;
        ide_spectrum.tilda_dSPECTRUM = tilda_dSPECTRUM;
        ide_spectrum.dSPECTRUM       = dSPECTRUM;
        ide_spectrum.dSPECTRUM_NO_MEAN = dSPECTRUM_NO_MEAN;
    end
    
%% Perform identification checks, i.e. find out which parameters are involved
    if checks_via_subsets
        % identification_checks_via_subsets is only for debugging
        [ide_dynamic, ide_reducedform, ide_moments, ide_spectrum, ide_minimal] = ...
            identification_checks_via_subsets(ide_dynamic, ide_reducedform, ide_moments, ide_spectrum, ide_minimal, totparam_nbr, modparam_nbr, options_ident, error_indicator);
         if ~error_indicator.identification_minimal
             ide_minimal.minimal_state_space=1;
         else
             ide_minimal.minimal_state_space=0;
         end
    else
        [ide_dynamic.cond, ide_dynamic.rank, ide_dynamic.ind0, ide_dynamic.indno, ide_dynamic.ino, ide_dynamic.Mco, ide_dynamic.Pco, ide_dynamic.jweak, ide_dynamic.jweak_pair] = ...
            identification_checks(dDYNAMIC(ind_dDYNAMIC,:)./norm_dDYNAMIC, 1, tol_rank, tol_sv, modparam_nbr);
        if ~options_ident.no_identification_reducedform && ~error_indicator.identification_reducedform
            [ide_reducedform.cond, ide_reducedform.rank, ide_reducedform.ind0, ide_reducedform.indno, ide_reducedform.ino, ide_reducedform.Mco, ide_reducedform.Pco, ide_reducedform.jweak, ide_reducedform.jweak_pair] = ...
                identification_checks(dREDUCEDFORM(ind_dREDUCEDFORM,:)./norm_dREDUCEDFORM, 1, tol_rank, tol_sv, totparam_nbr);
        end
        if ~options_ident.no_identification_moments && ~error_indicator.identification_moments
            [ide_moments.cond, ide_moments.rank, ide_moments.ind0, ide_moments.indno, ide_moments.ino, ide_moments.Mco, ide_moments.Pco, ide_moments.jweak, ide_moments.jweak_pair] = ...
                identification_checks(dMOMENTS(ind_dMOMENTS,:)./norm_dMOMENTS, 1, tol_rank, tol_sv, totparam_nbr);
        end
        if ~options_ident.no_identification_minimal 
            if ~error_indicator.identification_minimal
                [ide_minimal.cond, ide_minimal.rank, ide_minimal.ind0, ide_minimal.indno, ide_minimal.ino, ide_minimal.Mco, ide_minimal.Pco, ide_minimal.jweak, ide_minimal.jweak_pair] = ...
                    identification_checks(dMINIMAL(ind_dMINIMAL,:)./norm_dMINIMAL, 2, tol_rank, tol_sv, totparam_nbr);
                ide_minimal.minimal_state_space=1;
            else
                ide_minimal.minimal_state_space=0;
            end
        end
        if ~options_ident.no_identification_spectrum && ~error_indicator.identification_spectrum
            [ide_spectrum.cond, ide_spectrum.rank, ide_spectrum.ind0, ide_spectrum.indno, ide_spectrum.ino, ide_spectrum.Mco, ide_spectrum.Pco, ide_spectrum.jweak, ide_spectrum.jweak_pair] = ...
                identification_checks(tilda_dSPECTRUM, 3, tol_rank, tol_sv, totparam_nbr);
        end
    end
end
