function DERIVS = get_perturbation_params_derivs(M_, options_, estim_params_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state, indpmodel, indpstderr, indpcorr, d2flag)
% DERIVS = get_perturbation_params_derivs(M_, options_, estim_params_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state, indpmodel, indpstderr, indpcorr, d2flag)
% previously getH.m in dynare 4.5
% -------------------------------------------------------------------------
% Computes derivatives (with respect to parameters) of
%   (1) steady-state (ys) and covariance matrix of shocks (Sigma_e)
%   (2) dynamic model jacobians (g1, g2, g3)
%   (3) perturbation solution matrices:
%       * order==1: ghx,ghu
%       * order==2: ghx,ghu,ghxx,ghxu,ghuu,ghs2
%       * order==3: ghx,ghu,ghxx,ghxu,ghuu,ghs2,ghxxx,ghxxu,ghxuu,ghuuu,ghxss,ghuss
% Note that the order in the parameter Jacobians is the following:
% (1) stderr parameters (indpstderr)
% (2) corr parameters (indpcorr)
% (3) model parameters (indpmodel)
%
% =========================================================================
% INPUTS
%   M_:                     [structure]             storing the model information
%   options_:               [structure]             storing the options
%   estim_params_:          [structure]             storing the estimation information
%   dr                      [structure]             Reduced form model.
%   endo_steady_state       [vector]                steady state value for endogenous variables
%   exo_steady_state        [vector]                steady state value for exogenous variables
%   exo_det_steady_state    [vector]                steady state value for exogenous deterministic variables                                    
%   indpmodel:              [modparam_nbr by 1]     index of selected (estimated) parameters in M_.params;
%                                                   corresponds to model parameters (no stderr and no corr) 
%                                                   in estimated_params block
%   indpstderr:             [stderrparam_nbr by 1]  index of selected (estimated) standard errors,
%                                                   i.e. for all exogenous variables where 'stderr' is given 
%                                                   in the estimated_params block
%   indpcorr:               [corrparam_nbr by 2]    matrix of selected (estimated) correlations,
%                                                   i.e. for all exogenous variables where 'corr' is given in 
%                                                   the estimated_params block
%   d2flag:                 [boolean]               flag to compute second-order parameter derivatives of steady state 
%                                                   and first-order Kalman transition matrices
% -------------------------------------------------------------------------
% OUTPUTS
% DERIVS: Structure with the following fields:
%   dYss:       [endo_nbr by modparam_nbr] in DR order
%               Jacobian (wrt model parameters only) of steady state, i.e. ys(order_var,:)
%   dSigma_e:   [exo_nbr by exo_nbr by totparam_nbr] in declaration order
%                Jacobian (wrt to all paramters) of covariance matrix of shocks, i.e. Sigma_e
%   dg1:        [endo_nbr by yy0ex0_nbr by modparam_nbr] in DR order
%               Parameter Jacobian of first derivative (wrt dynamic model variables) of dynamic model (wrt to model parameters only)
%   dg2:        [endo_nbr by yy0ex0_nbr^2*modparam_nbr] in DR order
%               Parameter Jacobian of second derivative (wrt dynamic model variables) of dynamic model (wrt to model parameters only)
%               Note that instead of tensors we use matrix notation with blocks: dg2 = [dg2_dp1 dg2_dp2 ...],
%               where dg2_dpj is [endo_nbr by yy0ex0_nbr^2] and represents the derivative of g2 wrt parameter pj
%   dg3:        [endo_nbr by yy0ex0_nbr^3*modparam_nbr] in DR order
%               Parameter Jacobian of third derivative (wrt dynamic model variables) of dynamic model (wrt to model parameters only)
%               Note that instead of tensors we use matrix notation with blocks: dg3 = [dg3_dp1 dg3_dp2 ...],
%               where dg3_dpj is [endo_nbr by yy0ex0_nbr^3] and represents the derivative of g3 wrt parameter pj
%   dghx:       [endo_nbr by nspred by totparam_nbr] in DR order
%               Jacobian (wrt to all parameters) of first-order perturbation solution matrix ghx
%   dghu:       [endo_nbr by exo_nbr by totparam_nbr] in DR order
%               Jacobian (wrt to all parameters) of first-order perturbation solution matrix ghu
%   dOm:        [endo_nbr by endo_nbr by totparam_nbr] in DR order
%               Jacobian (wrt to all paramters) of Om = ghu*Sigma_e*transpose(ghu)
%   dghxx       [endo_nbr by nspred*nspred by totparam_nbr] in DR order
%               Jacobian (wrt to all parameters) of second-order perturbation solution matrix ghxx
%   dghxu       [endo_nbr by nspred*exo_nbr by totparam_nbr] in DR order
%               Jacobian (wrt to all parameters) of second-order perturbation solution matrix ghxu
%   dghuu       [endo_nbr by exo_nbr*exo_nbr by totparam_nbr] in DR order
%               Jacobian (wrt to all parameters) of second-order perturbation solution matrix ghuu
%   dghs2       [endo_nbr by totparam_nbr] in DR order
%               Jacobian (wrt to all parameters) of second-order perturbation solution matrix ghs2
%   dghxxx      [endo_nbr by nspred*nspred*nspred by totparam_nbr] in DR order
%               Jacobian (wrt to all parameters) of third-order perturbation solution matrix ghxxx
%   dghxxu      [endo_nbr by nspred*nspred*exo_nbr by totparam_nbr] in DR order
%               Jacobian (wrt to all parameters) of third-order perturbation solution matrix ghxxu
%   dghxuu      [endo_nbr by nspred*exo_nbr*exo_nbr by totparam_nbr] in DR order
%               Jacobian (wrt to all parameters) of third-order perturbation solution matrix ghxuu
%   dghuuu      [endo_nbr by exo_nbr*exo_nbr*exo_nbr by totparam_nbr] in DR order
%               Jacobian (wrt to all parameters) of third-order perturbation solution matrix ghuuu
%   dghxss      [endo_nbr by nspred by totparam_nbr] in DR order
%               Jacobian (wrt to all parameters) of third-order perturbation solution matrix ghxss
%   dghuss      [endo_nbr by exo_nbr by totparam_nbr] in DR order
%               Jacobian (wrt to all parameters) of third-order perturbation solution matrix ghuss
% if d2flag==true, we additional output:
%   d2KalmanA:  [endo_nbr*endo_nbr by totparam_nbr*(totparam_nbr+1)/2] in DR order
%               Unique entries of Hessian (wrt all parameters) of Kalman transition matrix A
%   d2Om:       [endo_nbr*(endo_nbr+1)/2 by totparam_nbr*(totparam_nbr+1)/2] in DR order
%               Unique entries of Hessian (wrt all parameters) of Om=ghu*Sigma_e*transpose(ghu)
%   d2Yss:      [endo_nbr by modparam_nbr by modparam_nbr] in DR order
%               Unique entries of Hessian (wrt model parameters only) of steady state ys(order_var,:)
%
% -------------------------------------------------------------------------
% This function is called by
%   * dsge_likelihood.m
%   * identification.get_jacobians.m
% -------------------------------------------------------------------------
% This function calls
%   * [fname,'.dynamic']
%   * [fname,'.dynamic_params_derivs']
%   * [fname,'.static']
%   * [fname,'.static_params_derivs']
%   * commutation
%   * dyn_vech
%   * dyn_unvech
%   * fjaco
%   * get_2nd_deriv (embedded)
%   * get_2nd_deriv_mat(embedded)
%   * get_all_parameters
%   * get_all_resid_2nd_derivs (embedded)
%   * get_hess_deriv (embedded)
%   * hessian_sparse
%   * sylvester3
%   * sylvester3a
%   * get_perturbation_params_derivs_numerical_objective
% =========================================================================
% Copyright © 2019-2020 Dynare Team
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
% Get fields from M_
Correlation_matrix = M_.Correlation_matrix;
dname              = M_.dname;
dynamic_tmp_nbr    = M_.dynamic_tmp_nbr;
endo_nbr           = M_.endo_nbr;
exo_nbr            = M_.exo_nbr;
exo_det_nbr        = M_.exo_det_nbr;
fname              = M_.fname;
lead_lag_incidence = M_.lead_lag_incidence;
nfwrd              = M_.nfwrd;
npred              = M_.npred;
nspred             = M_.nspred;
nstatic            = M_.nstatic;
params             = M_.params;
param_nbr          = M_.param_nbr;
Sigma_e            = M_.Sigma_e;

% Get fields from options_
analytic_derivation_mode = options_.analytic_derivation_mode;
    % analytic_derivation_mode: select method to compute Jacobians, default is 0
    % *  0: efficient sylvester equation method to compute analytical derivatives as in Ratto & Iskrev (2012)
    % *  1: kronecker products method to compute analytical derivatives as in Iskrev (2010), only for order=1
    % * -1: numerical two-sided finite difference method to compute numerical derivatives of all output arguments using function get_perturbation_params_derivs_numerical_objective.m
    % * -2: numerical two-sided finite difference method to compute numerically dYss, dg1, dg2, dg3, d2Yss and d2g1, the other output arguments are computed analytically as in kronflag=0
gstep        = options_.gstep;
order        = options_.order;
if isempty(options_.qz_criterium)
    % set default value for qz_criterium: if there are no unit roots one can use 1.0
    % If they are possible, you may have have multiple unit roots and the accuracy 
    % decreases when computing the eigenvalues in lyapunov_symm. Hence, we normally use 1+1e-6
    options_ = select_qz_criterium_value(options_);
end
qz_criterium = options_.qz_criterium;
threads_BC   = options_.threads.kronecker.sparse_hessian_times_B_kronecker_C;

% Get fields from dr
ghx = dr.ghx;
ghu = dr.ghu;
if order > 1
    ghxx = dr.ghxx;
    ghxu = dr.ghxu;
    ghuu = dr.ghuu;
    ghs2 = dr.ghs2;
end
if order > 2
    ghxxx = dr.ghxxx;
    ghxxu = dr.ghxxu;
    ghxuu = dr.ghxuu;
    ghuuu = dr.ghuuu;
    ghxss = dr.ghxss;
    ghuss = dr.ghuss;
end
order_var = dr.order_var;
ys        = dr.ys;

% Some checks
if exo_det_nbr > 0
    error('''get_perturbation_params_derivs'': not compatible with deterministic exogenous variables, please declare as endogenous.')
end
if order > 1 && analytic_derivation_mode == 1
    %analytic derivatives using Kronecker products is implemented only at first-order, at higher order we reset to analytic derivatives with sylvester equations
    %options_.analytic_derivation_mode = 0; fprintf('As order > 1, reset ''analytic_derivation_mode'' to 0\n');
    analytic_derivation_mode = 0; fprintf('As order > 1, reset ''analytic_derivation_mode'' to 0\n');
end

numerical_objective_fname = str2func('identification.get_perturbation_params_derivs_numerical_objective');
idx_states      = nstatic+(1:nspred); %index for state variables, in DR order
modparam_nbr    = length(indpmodel);  %number of selected model parameters
stderrparam_nbr = length(indpstderr); %number of selected stderr parameters
corrparam_nbr   = size(indpcorr,1);   %number of selected corr parameters
totparam_nbr    = modparam_nbr + stderrparam_nbr + corrparam_nbr; %total number of selected parameters
[I,~]           = find(lead_lag_incidence');                      %I is used to select nonzero columns of the Jacobian of endogenous variables in dynamic model files
yy0_nbr         = length(ys(I));                                  %number of dynamic variables
yy0ex0_nbr      = yy0_nbr+exo_nbr;                                %number of dynamic variables + exogenous variables
kyy0            = nonzeros(lead_lag_incidence(:,order_var)');     %index for nonzero entries in dynamic files at t-1,t,t+1 in DR order
kyy0ex0         = [kyy0; length(kyy0)+(1:exo_nbr)'];              %dynamic files include derivatives wrt exogenous variables, note that exo_det is always 0
if order > 1
    k2yy0ex0    = transpose(reshape(1:yy0ex0_nbr^2,yy0ex0_nbr,yy0ex0_nbr)); %index for the second dynamic derivatives, i.e. to evaluate the derivative of f wrt to yy0ex0(i) and yy0ex0(j), in DR order
end
if order > 2
    k3yy0ex0    = permute(reshape(transpose(reshape(1:yy0ex0_nbr^3,yy0ex0_nbr,yy0ex0_nbr^2)),yy0ex0_nbr,yy0ex0_nbr,yy0ex0_nbr),[2 1 3]); %index for the third dynamic derivative, i.e. df/(dyyex0_i*dyyex0_j*dyyex0_k)
end

% Check for purely backward or forward looking models
if size(lead_lag_incidence,1)<3
    if nfwrd == 0 %purely backward models
        klag       = lead_lag_incidence(1,order_var);        %indices of lagged (i.e. t-1) variables in dynamic files, columns are in DR order
        kcurr      = lead_lag_incidence(2,order_var);        %indices of current (i.e. t) variables in dynamic files, columns are in DR order
        klead      = zeros(1,size(lead_lag_incidence,2));          %indices of lead (i.e. t+1) variables in dynamic files, columns are in DR order
    elseif npred == 0 %purely forward models
        klag       = zeros(1,size(lead_lag_incidence,2));          %indices of lagged (i.e. t-1) variables in dynamic files, columns are in DR order
        kcurr      = lead_lag_incidence(1,order_var);        %indices of current (i.e. t) variables in dynamic files, columns are in DR order
        klead      = lead_lag_incidence(2,order_var);        %indices of lead (i.e. t+1) variables in dynamic files, columns are in DR order
    end
else %normal models
    klag       = lead_lag_incidence(1,order_var);        %indices of lagged (i.e. t-1) variables in dynamic files, columns are in DR order
    kcurr      = lead_lag_incidence(2,order_var);        %indices of current (i.e. t) variables in dynamic files, columns are in DR order
    klead      = lead_lag_incidence(3,order_var);        %indices of lead (i.e. t+1) variables in dynamic files, columns are in DR order
end

if analytic_derivation_mode < 0
    %Create auxiliary estim_params_ blocks if not available for numerical derivatives, estim_params_model contains only model parameters
    estim_params_model.np = length(indpmodel);
    estim_params_model.param_vals(:,1) = indpmodel;
    estim_params_model.nvx = 0; estim_params_model.ncx = 0; estim_params_model.nvn = 0; estim_params_model.ncn = 0;
    modparam1 = get_all_parameters(estim_params_model, M_);  %get all selected model parameters
    if ~isempty(indpstderr) && isempty(estim_params_.var_exo) %if there are stderr parameters but no estimated_params_block
        %provide temporary necessary information for stderr parameters
        estim_params_.nvx = length(indpstderr);
        estim_params_.var_exo(:,1) = indpstderr;
    end
    if ~isempty(indpcorr) && isempty(estim_params_.corrx) %if there are corr parameters but no estimated_params_block
        %provide temporary necessary information for stderr parameters
        estim_params_.ncx = size(indpcorr,1);
        estim_params_.corrx(:,1:2) = indpcorr;
    end
    if ~isfield(estim_params_,'nvn') %stderr of measurement errors not yet
        estim_params_.nvn = 0;
        estim_params_.var_endo = [];
    end
    if ~isfield(estim_params_,'ncn') %corr of measurement errors not yet
        estim_params_.ncn = 0;
        estim_params_.corrn = [];
    end
    if ~isempty(indpmodel) && isempty(estim_params_.param_vals) %if there are model parameters but no estimated_params_block
        %provide temporary necessary information for model parameters
        estim_params_.np = length(indpmodel);
        estim_params_.param_vals(:,1) = indpmodel;
    end
    xparam1 = get_all_parameters(estim_params_, M_);  %get all selected stderr, corr, and model parameters in estimated_params block, stderr and corr come first, then model parameters
end
if d2flag
    modparam_nbr2 = modparam_nbr*(modparam_nbr+1)/2; %number of unique entries of selected model parameters only in second-order derivative matrix
    totparam_nbr2 = totparam_nbr*(totparam_nbr+1)/2; %number of unique entries of all selected parameters in second-order derivative matrix
    %get indices of elements in second derivatives of parameters
    indp2tottot       = reshape(1:totparam_nbr^2,totparam_nbr,totparam_nbr);
    indp2stderrstderr = indp2tottot(1:stderrparam_nbr , 1:stderrparam_nbr);
    indp2stderrcorr   = indp2tottot(1:stderrparam_nbr , stderrparam_nbr+1:stderrparam_nbr+corrparam_nbr);
    indp2modmod       = indp2tottot(stderrparam_nbr+corrparam_nbr+1:stderrparam_nbr+corrparam_nbr+modparam_nbr , stderrparam_nbr+corrparam_nbr+1:stderrparam_nbr+corrparam_nbr+modparam_nbr);
    if totparam_nbr ~=1
        indp2tottot2 = dyn_vech(indp2tottot); %index of unique second-order derivatives
    else
        indp2tottot2 = indp2tottot;
    end
    if modparam_nbr ~= 1
        indp2modmod2 = dyn_vech(indp2modmod); %get rid of cross derivatives
    else
        indp2modmod2 = indp2modmod;
    end
    %Kalman transition matrices, as in kalman_transition_matrix.m
    KalmanA = zeros(endo_nbr,endo_nbr);
    KalmanA(:,idx_states) = ghx;
    KalmanB = ghu;
end

% Store some objects
DERIVS.indpmodel  = indpmodel;
DERIVS.indpstderr = indpstderr;
DERIVS.indpcorr   = indpcorr;

if analytic_derivation_mode == -1
%% numerical two-sided finite difference method using function get_perturbation_params_derivs_numerical_objective.m (previously thet2tau.m in Dynare 4.5) for
% Jacobian (wrt selected stderr, corr and model parameters) of
% - dynamic model derivatives: dg1, dg2, dg3
% - steady state (in DR order): dYss
% - covariance matrix of shocks: dSigma_e
% - perturbation solution matrices: dghx, dghu, dghxx, dghxu, dghuu, dghs2, dghxxx, dghxxu, dghxuu, dghuuu, dghxss, dghuss

    %Parameter Jacobian of covariance matrix and solution matrices (wrt selected stderr, corr and model paramters)
    dSig_gh         = identification.fjaco(numerical_objective_fname, xparam1, 'perturbation_solution', estim_params_, M_, options_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
    ind_Sigma_e     = (1:exo_nbr^2);
    ind_ghx         = ind_Sigma_e(end) + (1:endo_nbr*nspred);
    ind_ghu         = ind_ghx(end) + (1:endo_nbr*exo_nbr);
    DERIVS.dSigma_e = reshape(dSig_gh(ind_Sigma_e,:),[exo_nbr exo_nbr totparam_nbr]); %in tensor notation, wrt selected parameters
    DERIVS.dghx     = reshape(dSig_gh(ind_ghx,:),[endo_nbr nspred totparam_nbr]);     %in tensor notation, wrt selected parameters
    DERIVS.dghu     = reshape(dSig_gh(ind_ghu,:),[endo_nbr exo_nbr totparam_nbr]);    %in tensor notation, wrt selected parameters
    if order > 1
        ind_ghxx     = ind_ghu(end)  + (1:endo_nbr*nspred^2);
        ind_ghxu     = ind_ghxx(end) + (1:endo_nbr*nspred*exo_nbr);
        ind_ghuu     = ind_ghxu(end) + (1:endo_nbr*exo_nbr*exo_nbr);
        ind_ghs2     = ind_ghuu(end) + (1:endo_nbr);
        DERIVS.dghxx = reshape(dSig_gh(ind_ghxx,:), [endo_nbr nspred^2 totparam_nbr]);        %in tensor notation, wrt selected parameters
        DERIVS.dghxu = reshape(dSig_gh(ind_ghxu,:), [endo_nbr nspred*exo_nbr totparam_nbr]);  %in tensor notation, wrt selected parameters
        DERIVS.dghuu = reshape(dSig_gh(ind_ghuu,:), [endo_nbr exo_nbr*exo_nbr totparam_nbr]); %in tensor notation, wrt selected parameters
        DERIVS.dghs2 = reshape(dSig_gh(ind_ghs2,:), [endo_nbr totparam_nbr]);                 %in tensor notation, wrt selected parameters
    end
    if order > 2
        ind_ghxxx     = ind_ghs2(end)  + (1:endo_nbr*nspred^3);
        ind_ghxxu     = ind_ghxxx(end) + (1:endo_nbr*nspred^2*exo_nbr);
        ind_ghxuu     = ind_ghxxu(end) + (1:endo_nbr*nspred*exo_nbr^2);
        ind_ghuuu     = ind_ghxuu(end) + (1:endo_nbr*exo_nbr^3);
        ind_ghxss     = ind_ghuuu(end) + (1:endo_nbr*nspred);
        ind_ghuss     = ind_ghxss(end) + (1:endo_nbr*exo_nbr);
        DERIVS.dghxxx = reshape(dSig_gh(ind_ghxxx,:), [endo_nbr nspred^3 totparam_nbr]);         %in tensor notation, wrt selected parameters
        DERIVS.dghxxu = reshape(dSig_gh(ind_ghxxu,:), [endo_nbr nspred^2*exo_nbr totparam_nbr]); %in tensor notation, wrt selected parameters
        DERIVS.dghxuu = reshape(dSig_gh(ind_ghxuu,:), [endo_nbr nspred*exo_nbr^2 totparam_nbr]); %in tensor notation, wrt selected parameters
        DERIVS.dghuuu = reshape(dSig_gh(ind_ghuuu,:), [endo_nbr exo_nbr^3 totparam_nbr]);        %in tensor notation, wrt selected parameters
        DERIVS.dghxss = reshape(dSig_gh(ind_ghxss,:), [endo_nbr nspred totparam_nbr]);           %in tensor notation, wrt selected parameters
        DERIVS.dghuss = reshape(dSig_gh(ind_ghuss,:), [endo_nbr exo_nbr totparam_nbr]);          %in tensor notation, wrt selected parameters
    end
    % Parameter Jacobian of Om=ghu*Sigma_e*ghu' and Correlation_matrix (wrt selected stderr, corr and model paramters)
    DERIVS.dOm                 = zeros(endo_nbr,endo_nbr,totparam_nbr); %initialize in tensor notation
    DERIVS.dCorrelation_matrix = zeros(exo_nbr,exo_nbr,totparam_nbr);   %initialize in tensor notation
    if ~isempty(indpstderr) %derivatives of ghu wrt stderr parameters are zero by construction
        for jp=1:stderrparam_nbr
            DERIVS.dOm(:,:,jp) = ghu*DERIVS.dSigma_e(:,:,jp)*ghu';
        end
    end
    if ~isempty(indpcorr)  %derivatives of ghu wrt corr parameters are zero by construction
        for jp=1:corrparam_nbr
            DERIVS.dOm(:,:,stderrparam_nbr+jp) = ghu*DERIVS.dSigma_e(:,:,stderrparam_nbr+jp)*ghu';
            DERIVS.dCorrelation_matrix(indpcorr(jp,1),indpcorr(jp,2),stderrparam_nbr+jp) = 1;
            DERIVS.dCorrelation_matrix(indpcorr(jp,2),indpcorr(jp,1),stderrparam_nbr+jp) = 1;%symmetry
        end
    end
    if ~isempty(indpmodel)  %derivatives of Sigma_e wrt model parameters are zero by construction
        for jp=1:modparam_nbr
            DERIVS.dOm(:,:,stderrparam_nbr+corrparam_nbr+jp) = DERIVS.dghu(:,:,stderrparam_nbr+corrparam_nbr+jp)*Sigma_e*ghu' + ghu*Sigma_e*DERIVS.dghu(:,:,stderrparam_nbr+corrparam_nbr+jp)';
        end
    end

    %Parameter Jacobian of dynamic model derivatives (wrt selected model parameters only)
    dYss_g = identification.fjaco(numerical_objective_fname, modparam1, 'dynamic_model', estim_params_model, M_, options_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
    ind_Yss = 1:endo_nbr;
    if options_.discretionary_policy || options_.ramsey_policy
        ind_g1 = ind_Yss(end) + (1:M_.eq_nbr*yy0ex0_nbr);
    else
        ind_g1 = ind_Yss(end) + (1:endo_nbr*yy0ex0_nbr);
    end
    DERIVS.dYss = dYss_g(ind_Yss, :); %in tensor notation, wrt selected model parameters only
    if options_.discretionary_policy || options_.ramsey_policy
        DERIVS.dg1 = reshape(dYss_g(ind_g1,:),[M_.eq_nbr, yy0ex0_nbr, modparam_nbr]); %in tensor notation, wrt selected model parameters only
    else
        DERIVS.dg1 = reshape(dYss_g(ind_g1,:),[endo_nbr, yy0ex0_nbr, modparam_nbr]); %in tensor notation, wrt selected model parameters only
    end
    if order > 1
        ind_g2 = ind_g1(end) + (1:endo_nbr*yy0ex0_nbr^2);
        DERIVS.dg2 = reshape(sparse(dYss_g(ind_g2,:)),[endo_nbr, yy0ex0_nbr^2*modparam_nbr]); %blockwise in matrix notation, i.e. [dg2_dp1 dg2_dp2 ...], where dg2_dpj has dimension endo_nbr by yy0ex0_nbr^2
    end
    if order > 2
        ind_g3 = ind_g2(end) + (1:endo_nbr*yy0ex0_nbr^3);
        DERIVS.dg3 = reshape(sparse(dYss_g(ind_g3,:)),[endo_nbr, yy0ex0_nbr^3*modparam_nbr]); %blockwise in matrix notation, i.e. [dg3_dp1 dg3_dp2 ...], where dg3_dpj has dimension endo_nbr by yy0ex0_nbr^3
    end

    if d2flag
        % Hessian (wrt paramters) of steady state and first-order solution matrices ghx and Om
        % note that hessian_sparse.m (contrary to hessian.m) does not take symmetry into account, but focuses already on unique values
        options_.order = 1; %make sure only first order
        d2Yss_KalmanA_Om = identification.hessian_sparse(numerical_objective_fname, xparam1, gstep, 'Kalman_Transition', estim_params_, M_, options_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
        options_.order = order; %make sure to set back
        ind_KalmanA = ind_Yss(end) + (1:endo_nbr^2);
        DERIVS.d2KalmanA = d2Yss_KalmanA_Om(ind_KalmanA, indp2tottot2);               %only unique elements
        DERIVS.d2Om      = d2Yss_KalmanA_Om(ind_KalmanA(end)+1:end , indp2tottot2);   %only unique elements
        DERIVS.d2Yss     = zeros(endo_nbr,modparam_nbr,modparam_nbr);                 %initialize
        for j = 1:endo_nbr
            DERIVS.d2Yss(j,:,:) = dyn_unvech(full(d2Yss_KalmanA_Om(j,indp2modmod2))); %Hessian for d2Yss, but without cross derivatives
        end
    end

    return %[END OF MAIN FUNCTION]!!!!!
end

if analytic_derivation_mode == -2
%% Numerical two-sided finite difference method to compute parameter derivatives of steady state and dynamic model,
% i.e. dYss, dg1, dg2, dg3 as well as d2Yss, d2g1 numerically.
% The parameter derivatives of perturbation solution matrices are computed analytically below (analytic_derivation_mode=0)
    if order == 3
        [~, g1, g2, g3] = feval([fname,'.dynamic'], ys(I), exo_steady_state', params, ys, 1);
        g3 = identification.unfold_g3(g3, yy0ex0_nbr);
    elseif order == 2
        [~, g1, g2] = feval([fname,'.dynamic'], ys(I), exo_steady_state', params, ys, 1);
    elseif order == 1
        [~, g1] = feval([fname,'.dynamic'], ys(I), exo_steady_state', params, ys, 1);
    end

    if d2flag
        % computation of d2Yss and d2g1
        % note that hessian_sparse does not take symmetry into account, i.e. compare hessian_sparse.m to hessian.m, but focuses already on unique values, which are duplicated below
        options_.order = 1; %d2flag requires only first order
        d2Yss_g1 = identification.hessian_sparse(numerical_objective_fname, modparam1, gstep, 'dynamic_model', estim_params_model, M_, options_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state);  % d2flag requires only first-order
        options_.order = order; %make sure to set back the order
        d2Yss = reshape(full(d2Yss_g1(1:endo_nbr,:)), [endo_nbr modparam_nbr modparam_nbr]); %put into tensor notation
        for j=1:endo_nbr
            d2Yss(j,:,:) = dyn_unvech(dyn_vech(d2Yss(j,:,:))); %add duplicate values to full hessian
        end
        d2g1_full = d2Yss_g1(endo_nbr+1:end,:);
        %store only nonzero unique entries and the corresponding indices of d2g1:
        %  rows: respective derivative term
        %  1st column: equation number of the term appearing
        %  2nd column: column number of variable in Jacobian of the dynamic model
        %  3rd column: number of the first parameter in derivative
        %  4th column: number of the second parameter in derivative
        %  5th column: value of the Hessian term
        ind_d2g1 = find(d2g1_full);
        d2g1 = zeros(length(ind_d2g1),5);
        for j=1:length(ind_d2g1)
            [i1, i2] = ind2sub(size(d2g1_full),ind_d2g1(j));
            [ig1, ig2] = ind2sub(size(g1),i1);
            [ip1, ip2] = ind2sub([modparam_nbr modparam_nbr],i2);
            d2g1(j,:) = [ig1 ig2 ip1 ip2 d2g1_full(ind_d2g1(j))];
        end
        clear d2g1_full d2Yss_g1;
    end

    %Parameter Jacobian of dynamic model derivatives (wrt selected model parameters only)
    dYss_g  = identification.fjaco(numerical_objective_fname, modparam1, 'dynamic_model', estim_params_model, M_, options_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
    ind_Yss = 1:endo_nbr;
    ind_g1  = ind_Yss(end) + (1:endo_nbr*yy0ex0_nbr);
    dYss    = dYss_g(ind_Yss,:); %in tensor notation, wrt selected model parameters only
    dg1     = reshape(dYss_g(ind_g1,:),[endo_nbr,yy0ex0_nbr,modparam_nbr]); %in tensor notation
    if order > 1
        ind_g2 = ind_g1(end) + (1:endo_nbr*yy0ex0_nbr^2);
        dg2 = reshape(sparse(dYss_g(ind_g2,:)),[endo_nbr, yy0ex0_nbr^2*modparam_nbr]); %blockwise in matrix notation, i.e. [dg2_dp1 dg2_dp2 ...], where dg2_dpj has dimension endo_nbr by yy0ex0_nbr^2
    end
    if order > 2
        ind_g3 = ind_g2(end) + (1:endo_nbr*yy0ex0_nbr^3);
        dg3 = reshape(sparse(dYss_g(ind_g3,:)), [endo_nbr, yy0ex0_nbr^3*modparam_nbr]); %blockwise in matrix notation, i.e. [dg3_dp1 dg3_dp2 ...], where dg3_dpj has dimension endo_nbr by yy0ex0_nbr^3
    end
    clear dYss_g

elseif (analytic_derivation_mode == 0 || analytic_derivation_mode == 1)
    if ~exist(['+' fname filesep 'static_params_derivs.m'],'file')
        error('For analytical parameter derivatives ''static_params_derivs.m'' file is needed, this can be created by putting identification(order=%d) into your mod file.',order)
    end
    if ~exist(['+' fname filesep 'dynamic_params_derivs.m'],'file')
        error('For analytical parameter derivatives ''dynamic_params_derivs.m'' file is needed, this can be created by putting identification(order=%d) into your mod file.',order)
    end
    %% Analytical computation of Jacobian and Hessian (wrt selected model parameters) of steady state, i.e. dYss and d2Yss
    [~, g1_static] = feval([fname,'.static'], ys, exo_steady_state', params); %g1_static is [endo_nbr by endo_nbr] first-derivative (wrt all endogenous variables) of static model equations f, i.e. df/dys, in declaration order
    rp_static = feval([fname,'.static_params_derivs'], ys, exo_steady_state', params); %rp_static is [endo_nbr by param_nbr] first-derivative (wrt all model parameters) of static model equations f, i.e. df/dparams, in declaration order
    dys = -g1_static\rp_static; %use implicit function theorem (equation 5 of Ratto and Iskrev (2012) to compute [endo_nbr by param_nbr] first-derivative (wrt all model parameters) of steady state for all endogenous variables analytically, note that dys is in declaration order
    d2ys = zeros(endo_nbr, param_nbr, param_nbr); %initialize in tensor notation, note that d2ys is only needed for d2flag, i.e. for g1pp
    if d2flag
        [~, ~, g2_static] = feval([fname,'.static'], ys, exo_steady_state', params); %g2_static is [endo_nbr by endo_nbr^2] second derivative (wrt all endogenous variables) of static model equations f, i.e. d(df/dys)/dys, in declaration order
        if order < 3
            [~, g1, g2, g3] = feval([fname,'.dynamic'], ys(I), exo_steady_state', params, ys, 1); %note that g3 does not contain symmetric elements
            g3 = identification.unfold_g3(g3, yy0ex0_nbr); %add symmetric elements to g3
        else
            T  = NaN(sum(dynamic_tmp_nbr(1:5)));
            T  = feval([fname, '.dynamic_g4_tt'], T, ys(I), exo_steady_state', params, ys, 1);
            g1 = feval([fname, '.dynamic_g1'],    T, ys(I), exo_steady_state', params, ys, 1, false); %g1 is [endo_nbr by yy0ex0_nbr first derivative (wrt all dynamic variables) of dynamic model equations, i.e. df/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
            g2 = feval([fname, '.dynamic_g2'],    T, ys(I), exo_steady_state', params, ys, 1, false); %g2 is [endo_nbr by yy0ex0_nbr^2] second derivative (wrt all dynamic variables) of dynamic model equations, i.e. d(df/dyy0ex0)/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
            g3 = feval([fname, '.dynamic_g3'],    T, ys(I), exo_steady_state', params, ys, 1, false); %note that g3 does not contain symmetric elements
            g4 = feval([fname, '.dynamic_g4'],    T, ys(I), exo_steady_state', params, ys, 1, false); %note that g4 does not contain symmetric elements
            g3 = identification.unfold_g3(g3, yy0ex0_nbr); %add symmetric elements to g3, %g3 is [endo_nbr by yy0ex0_nbr^3] third-derivative (wrt all dynamic variables) of dynamic model equations, i.e. (d(df/dyy0ex0)/dyy0ex0)/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
            g4 = identification.unfold_g4(g4, yy0ex0_nbr); %add symmetric elements to g4, %g4 is [endo_nbr by yy0ex0_nbr^4] fourth-derivative (wrt all dynamic variables) of dynamic model equations, i.e. ((d(df/dyy0ex0)/dyy0ex0)/dyy0ex0)/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
        end
        %g1 is [endo_nbr by yy0ex0_nbr first derivative (wrt all dynamic variables) of dynamic model equations, i.e. df/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
        %g2 is [endo_nbr by yy0ex0_nbr^2] second derivative (wrt all dynamic variables) of dynamic model equations, i.e. d(df/dyy0ex0)/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
        %g3 is [endo_nbr by yy0ex0_nbr^3] third-derivative (wrt all dynamic variables) of dynamic model equations, i.e. (d(df/dyy0ex0)/dyy0ex0)/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
        [~, g1p_static, rpp_static] = feval([fname,'.static_params_derivs'], ys, exo_steady_state', params);
        %g1p_static is [endo_nbr by endo_nbr by param_nbr] first derivative (wrt all model parameters) of first-derivative (wrt all endogenous variables) of static model equations f, i.e. (df/dys)/dparams, in declaration order
        %rpp_static is [#second_order_residual_terms by 4] and contains nonzero values and corresponding indices of second derivatives (wrt all model parameters) of static model equations f, i.e. d(df/dparams)/dparams, in declaration order, where
        %           column 1 contains equation number; column 2 contains first parameter; column 3 contains second parameter; column 4 contains value of derivative
        rpp_static = get_all_resid_2nd_derivs(rpp_static, endo_nbr, param_nbr); %make full matrix out of nonzero values and corresponding indices
            %rpp_static is [endo_nbr by param_nbr by param_nbr] second derivatives (wrt all model parameters) of static model equations, i.e. d(df/dparams)/dparams, in declaration order
        if isempty(find(g2_static))
            %auxiliary expression on page 8 of Ratto and Iskrev (2012) is zero, i.e. gam = 0
            for j = 1:param_nbr
                %using the implicit function theorem, equation 15 on page 7 of Ratto and Iskrev (2012)
                d2ys(:,:,j) = -g1_static\rpp_static(:,:,j);
                %d2ys is [endo_nbr by param_nbr by param_nbr] second-derivative (wrt all model parameters) of steady state of all endogenous variables, i.e. d(dys/dparams)/dparams, in declaration order
            end
        else
            gam = zeros(endo_nbr,param_nbr,param_nbr); %initialize auxiliary expression on page 8 of Ratto and Iskrev (2012)
            for j = 1:endo_nbr
                tmp_g1p_static_dys = (squeeze(g1p_static(j,:,:))'*dys);
                gam(j,:,:) = transpose(reshape(g2_static(j,:),[endo_nbr endo_nbr])*dys)*dys + tmp_g1p_static_dys + tmp_g1p_static_dys';
            end
            for j = 1:param_nbr
                %using the implicit function theorem, equation 15 on page 7 of Ratto and Iskrev (2012)
                d2ys(:,:,j) = -g1_static\(rpp_static(:,:,j)+gam(:,:,j));
                %d2ys is [endo_nbr by param_nbr by param_nbr] second-derivative (wrt all model parameters) of steady state of all endogenous variables, i.e. d(dys/dparams)/dparams, in declaration order
            end
            clear g1p_static g2_static tmp_g1p_static_dys gam
        end
    end
    %handling of steady state for nonstationary variables
    if any(any(isnan(dys)))
        [U,T] = schur(g1_static);
        e1 = abs(ordeig(T)) < qz_criterium-1;
        k = sum(e1);       % Number of non stationary variables.
                           % Number of stationary variables: n = length(e1)-k
        [U,T] = ordschur(U,T,e1);
        T = T(k+1:end,k+1:end);
        %using implicit function theorem, equation 5 of Ratto and Iskrev (2012), in declaration order
        dys = -U(:,k+1:end)*(T\U(:,k+1:end)')*rp_static;
        if d2flag
            fprintf('Computation of d2ys for nonstationary variables is not yet correctly handled if g2_static is nonempty, but continue anyways...\n')
            for j = 1:param_nbr
                %using implicit function theorem, equation 15 of Ratto and Iskrev (2012), in declaration order
                d2ys(:,:,j) = -U(:,k+1:end)*(T\U(:,k+1:end)')*rpp_static(:,:,j); %THIS IS NOT CORRECT, IF g2_static IS NONEMPTY. WE NEED TO ADD GAM [willi]
            end
        end
    end

    if d2flag
        if order < 3
            [~, g1p, ~, g1pp, g2p] = feval([fname,'.dynamic_params_derivs'], ys(I), exo_steady_state', params, ys, 1, dys, d2ys);
        else
            [~, g1p, ~, g1pp, g2p, g3p] = feval([fname,'.dynamic_params_derivs'], ys(I), exo_steady_state', params, ys, 1, dys, d2ys);
        end
        %g1pp are nonzero values and corresponding indices of second-derivatives (wrt all model parameters) of first-derivative (wrt all dynamic variables) of dynamic model equations, i.e. d(d(df/dyy0ex0)/dparam)/dparam, rows are in declaration order, first column in declaration order
        d2Yss = d2ys(order_var,indpmodel,indpmodel); %[endo_nbr by mod_param_nbr by mod_param_nbr], put into DR order and focus only on selected model parameters
    else
        if order == 1
            [~, g1p] = feval([fname,'.dynamic_params_derivs'], ys(I), exo_steady_state', params, ys, 1, dys, d2ys);
            %g1p is [endo_nbr by yy0ex0_nbr by param_nbr] first-derivative (wrt all model parameters) of first-derivative (wrt all dynamic variables) of dynamic model equations, i.e. d(df/dyy0ex0)/dparam, rows are in declaration order, column in lead_lag_incidence order
            [~, g1, g2 ] = feval([fname,'.dynamic'], ys(I), exo_steady_state', params, ys, 1);
                %g1 is [endo_nbr by yy0ex0_nbr first derivative (wrt all dynamic variables) of dynamic model equations, i.e. df/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
                %g2 is [endo_nbr by yy0ex0_nbr^2] second derivatives (wrt all dynamic variables) of dynamic model equations, i.e. d(df/dyy0ex0)/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
        elseif order == 2
            [~, g1p, ~, ~, g2p] = feval([fname,'.dynamic_params_derivs'], ys(I), exo_steady_state', params, ys, 1, dys, d2ys);
            %g1p is [endo_nbr by yy0ex0_nbr by param_nbr] first-derivative (wrt all model parameters) of first-derivative (wrt all dynamic variables) of dynamic model equations, i.e. d(df/dyy0ex0)/dparam, rows are in declaration order, column in lead_lag_incidence order
            %g2p are nonzero values and corresponding indices of first-derivative (wrt all model parameters) of second-derivatives (wrt all dynamic variables) of dynamic model equations, i.e. d(d(df/dyy0ex0)/dyy0ex0)/dparam, rows are in declaration order, first and second column in declaration order
            [~, g1, g2, g3] = feval([fname,'.dynamic'], ys(I), exo_steady_state', params, ys, 1); %note that g3 does not contain symmetric elements
            g3 = identification.unfold_g3(g3, yy0ex0_nbr); %add symmetric elements to g3
                %g1 is [endo_nbr by yy0ex0_nbr first derivative (wrt all dynamic variables) of dynamic model equations, i.e. df/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
                %g2 is [endo_nbr by yy0ex0_nbr^2] second derivative (wrt all dynamic variables) of dynamic model equations, i.e. d(df/dyy0ex0)/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
                %g3 is [endo_nbr by yy0ex0_nbr^3] third-derivative (wrt all dynamic variables) of dynamic model equations, i.e. (d(df/dyy0ex0)/dyy0ex0)/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
        elseif order == 3
            [~, g1p, ~, ~, g2p, g3p] = feval([fname,'.dynamic_params_derivs'], ys(I), exo_steady_state', params, ys, 1, dys, d2ys);
            %g1p is [endo_nbr by yy0ex0_nbr by param_nbr] first-derivative (wrt all model parameters) of first-derivative (wrt all dynamic variables) of dynamic model equations, i.e. d(df/dyy0ex0)/dparam, rows are in declaration order, column in lead_lag_incidence order
            %g2p are nonzero values and corresponding indices of first-derivative (wrt all model parameters) of second-derivatives (wrt all dynamic variables) of dynamic model equations, i.e. d(d(df/dyy0ex0)/dyy0ex0)/dparam, rows are in declaration order, first and second column in declaration order
            %g3p are nonzero values and corresponding indices of first-derivative (wrt all model parameters) of third-derivatives (wrt all dynamic variables) of dynamic model equations, i.e. d(d(d(df/dyy0ex0)/dyy0ex0)/dyy0ex0)/dparam, rows are in declaration order, first, second and third column in declaration order
            T  = NaN(sum(dynamic_tmp_nbr(1:5)));
            T  = feval([fname, '.dynamic_g4_tt'], T, ys(I), exo_steady_state', params, ys, 1);
            g1 = feval([fname, '.dynamic_g1'],    T, ys(I), exo_steady_state', params, ys, 1, false); %g1 is [endo_nbr by yy0ex0_nbr first derivative (wrt all dynamic variables) of dynamic model equations, i.e. df/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
            g2 = feval([fname, '.dynamic_g2'],    T, ys(I), exo_steady_state', params, ys, 1, false); %g2 is [endo_nbr by yy0ex0_nbr^2] second derivative (wrt all dynamic variables) of dynamic model equations, i.e. d(df/dyy0ex0)/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
            g3 = feval([fname, '.dynamic_g3'],    T, ys(I), exo_steady_state', params, ys, 1, false); %note that g3 does not contain symmetric elements
            g4 = feval([fname, '.dynamic_g4'],    T, ys(I), exo_steady_state', params, ys, 1, false); %note that g4 does not contain symmetric elements
            g3 = identification.unfold_g3(g3, yy0ex0_nbr); %add symmetric elements to g3, %g3 is [endo_nbr by yy0ex0_nbr^3] third-derivative (wrt all dynamic variables) of dynamic model equations, i.e. (d(df/dyy0ex0)/dyy0ex0)/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
            g4 = identification.unfold_g4(g4, yy0ex0_nbr); %add symmetric elements to g4, %g4 is [endo_nbr by yy0ex0_nbr^4] fourth-derivative (wrt all dynamic variables) of dynamic model equations, i.e. ((d(df/dyy0ex0)/dyy0ex0)/dyy0ex0)/dyy0ex0, rows are in declaration order, columns in lead_lag_incidence order
        end
    end
    % Parameter Jacobian of steady state in different orderings, note dys is in declaration order
    dYss = dys(order_var, indpmodel); %in DR-order, focus only on selected model parameters
    dyy0 = dys(I,:); %in lead_lag_incidence order, Jacobian of dynamic (without exogenous) variables, focus on all model parameters
    dyy0ex0 = sparse([dyy0; zeros(exo_nbr,param_nbr)]); %in lead_lag_incidence order, Jacobian of dynamic (with exogenous) variables, focus on all model parameters

    %% Analytical computation of Jacobian (wrt selected model parameters) of first-derivative of dynamic model, i.e. dg1
    %  Note that we use the implicit function theorem (see Ratto and Iskrev (2012) formula 7):
    %   Let g1 = df/dyy0ex0 be the first derivative (wrt all dynamic variables) of the dynamic model, then
    %   dg1 denotes the first-derivative (wrt model parameters) of g1 evaluated at the steady state.
    %   Note that g1 is a function of both the model parameters p and of the steady state of all dynamic variables, which also depend on the model parameters.
    %   Hence, implicitly g1=g1(p,yy0ex0(p)) and dg1 consists of two parts:
    %   (1) g1p: direct derivative (wrt to all model parameters) given by the preprocessor
    %   (2) g2*dyy0ex0: contribution of derivative of steady state of dynamic variables (wrt all model parameters)
    %   Note that in a stochastic context ex0 and dex0 is always zero and hence can be skipped in the computations
    dg1 = zeros(endo_nbr,yy0ex0_nbr,param_nbr); %initialize part (2)
    for j = 1:endo_nbr
        [II, JJ] = ind2sub([yy0ex0_nbr yy0ex0_nbr], find(g2(j,:))); %g2 is [endo_nbr by yy0ex0_nbr^2]
        for i = 1:yy0ex0_nbr
            is = find(II==i);
            is = is(find(JJ(is)<=yy0_nbr)); %focus only on dr.ys(I) derivatives as exogenous variables are 0 in a stochastic context
            if ~isempty(is)
                tmp_g2 = full(g2(j,find(g2(j,:))));
                dg1(j,i,:) = tmp_g2(is)*dyy0(JJ(is),:); %put into tensor notation
            end
        end
    end
    dg1 = g1p + dg1; %add part (1)
    dg1 = dg1(:,:,indpmodel); %focus only on selected model parameters

    if order>1
        %% Analytical computation of Jacobian (wrt selected model parameters) of second derivative of dynamic model, i.e. dg2
        %  We use the implicit function theorem:
        %   Let g2 = d2f/(dyy0ex0*dyy0ex0) denote the second derivative (wrt all dynamic variables) of the dynamic model, then
        %   dg2 denotes the first-derivative (wrt all model parameters) of g2 evaluated at the steady state.
        %   Note that g2 is a function of both the model parameters p and of the steady state of all dynamic variables, which also depend on the parameters.
        %   Hence, implicitly g2=g2(p,yy0ex0(p)) and dg2 consists of two parts:
        %   (1) g2p: direct derivative wrt to all model parameters given by the preprocessor
        %   and
        %   (2) g3*dyy0ex0: contribution of derivative of steady state of dynamic variables (wrt all model parameters)
        % Note that we focus on selected model parameters only, i.e. indpmodel
        % Also note that we stack the parameter derivatives blockwise instead of in tensors
        dg2 = reshape(g3,[endo_nbr*yy0ex0_nbr^2 yy0ex0_nbr])*dyy0ex0(:,indpmodel); %part (2)
        dg2 = reshape(dg2, [endo_nbr, yy0ex0_nbr^2*modparam_nbr]);
        for jj = 1:size(g2p,1)
            jpos = find(indpmodel==g2p(jj,4));
            if jpos~=0
                dg2(g2p(jj,1), k2yy0ex0(g2p(jj,2),g2p(jj,3))+(jpos-1)*yy0ex0_nbr^2) = dg2(g2p(jj,1), k2yy0ex0(g2p(jj,2),g2p(jj,3))+(jpos-1)*yy0ex0_nbr^2) + g2p(jj,5); %add part (1)
            end
        end
    end

    if order>2
        %% Analytical computation of Jacobian (wrt selected model parameters) of third derivative of dynamic model, i.e. dg3
        %  We use the implicit function theorem:
        %   Let g3 = d3f/(dyy0ex0*dyy0ex0*dyy0ex0) denote the third derivative (wrt all dynamic variables) of the dynamic model, then
        %   dg3 denotes the first-derivative (wrt all model parameters) of g3 evaluated at the steady state.
        %   Note that g3 is a function of both the model parameters p and of the steady state of all dynamic variables, which also depend on the parameters.
        %   Hence, implicitly g3=g3(p,yy0ex0(p)) and dg3 consists of two parts:
        %   (1) g3p: direct derivative wrt to all model parameters given by the preprocessor
        %   and
        %   (2) g4*dyy0ex0: contribution of derivative of steady state of dynamic variables (wrt all model parameters)
        % Note that we focus on selected model parameters only, i.e. indpmodel
        % Also note that we stack the parameter derivatives blockwise instead of in tensors
        dg3 = reshape(g4,[endo_nbr*yy0ex0_nbr^3 yy0ex0_nbr])*dyy0ex0(:,indpmodel); %part (2)
        dg3 = reshape(dg3, [endo_nbr, yy0ex0_nbr^3*modparam_nbr]);
        for jj = 1:size(g3p,1)
            jpos = find(indpmodel==g3p(jj,5));
            if jpos~=0
                idyyy = unique(perms([g3p(jj,2) g3p(jj,3) g3p(jj,4)]),'rows'); %note that g3p does not contain symmetric terms, so we use the perms and unique functions
                for k = 1:size(idyyy,1)
                    dg3(g3p(jj,1), k3yy0ex0(idyyy(k,1),idyyy(k,2),idyyy(k,3))+(jpos-1)*yy0ex0_nbr^3) = dg3(g3p(jj,1), k3yy0ex0(idyyy(k,1),idyyy(k,2),idyyy(k,3))+(jpos-1)*yy0ex0_nbr^3) + g3p(jj,6); %add part (1)
                end
            end
        end
    end

    if d2flag
        %% Analytical computation of Hessian (wrt selected model parameters) of first-derivative of dynamic model, i.e. d2g1
        % We use the implicit function theorem (above we already computed: dg1 = g1p + g2*dyy0ex0):
        % Accordingly, d2g1, the second-derivative (wrt model parameters), consists of five parts (ignoring transposes, see Ratto and Iskrev (2012) formula 16)
        %   (1) d(g1p)/dp                             = g1pp
        %   (2) d(g1p)/dyy0ex0*d(yy0ex0)/dp           = g2p * dyy0ex0
        %   (3) d(g2)/dp * dyy0ex0                    = g2p * dyy0ex0
        %   (4) d(g2)/dyy0ex0*d(dyy0ex0)/dp * dyy0ex0 = g3  * dyy0ex0 * dyy0ex0
        %   (5) g2 * d(dyy0ex0)/dp                    = g2  * d2yy0ex0
        %   Note that part 2 and 3 are equivalent besides the use of transpose (see Ratto and Iskrev (2012) formula 16)
        d2g1_full = sparse(endo_nbr*yy0ex0_nbr, param_nbr*param_nbr);  %initialize
        g3_tmp = reshape(g3,[endo_nbr*yy0ex0_nbr*yy0ex0_nbr yy0ex0_nbr]);
        d2g1_part4_left = sparse(endo_nbr*yy0ex0_nbr*yy0ex0_nbr,param_nbr);
        for j = 1:param_nbr
            %compute first two terms of part 4
            d2g1_part4_left(:,j) = g3_tmp*dyy0ex0(:,j);
        end
        for j=1:endo_nbr
            %Note that in the following we skip exogenous variables as these are 0 by construction in a stochastic setting
            d2g1_part5 = reshape(g2(j,:), [yy0ex0_nbr yy0ex0_nbr]);
            d2g1_part5 = d2g1_part5(:,1:yy0_nbr)*reshape(d2ys(I,:,:),[yy0_nbr,param_nbr*param_nbr]);
            for i=1:yy0ex0_nbr
                ind_part4 = sub2ind([endo_nbr yy0ex0_nbr yy0ex0_nbr], ones(yy0ex0_nbr,1)*j ,ones(yy0ex0_nbr,1)*i, (1:yy0ex0_nbr)');
                d2g1_part4 = (d2g1_part4_left(ind_part4,:))'*dyy0ex0;
                d2g1_part2_and_part3 = (get_hess_deriv(g2p,j,i,yy0ex0_nbr,param_nbr))'*dyy0ex0;
                d2g1_part1 = get_2nd_deriv_mat(g1pp,j,i,param_nbr);
                d2g1_tmp = d2g1_part1 + d2g1_part2_and_part3 + d2g1_part2_and_part3' + d2g1_part4 + reshape(d2g1_part5(i,:,:),[param_nbr param_nbr]);
                d2g1_tmp = d2g1_tmp(indpmodel,indpmodel); %focus only on model parameters
                if any(any(d2g1_tmp))
                    ind_d2g1_tmp = find(triu(d2g1_tmp));
                    d2g1_full(sub2ind([endo_nbr yy0ex0_nbr],j,i), ind_d2g1_tmp) = transpose(d2g1_tmp(ind_d2g1_tmp));
                end
            end
        end
        clear d2g1_tmp d2g1_part1 d2g1_part2_and_part3 d2g1_part4 d2g1_part4_left d2g1_part5
        %store only nonzero entries and the corresponding indices of d2g1:
        %  1st column: equation number of the term appearing
        %  2nd column: column number of variable in Jacobian of the dynamic model
        %  3rd column: number of the first parameter in derivative
        %  4th column: number of the second parameter in derivative
        %  5th column: value of the Hessian term
        ind_d2g1 = find(d2g1_full);
        d2g1 = zeros(length(ind_d2g1),5);
        for j=1:length(ind_d2g1)
            [i1, i2] = ind2sub(size(d2g1_full),ind_d2g1(j));
            [ig1, ig2] = ind2sub(size(g1),i1);
            [ip1, ip2] = ind2sub([modparam_nbr modparam_nbr],i2);
            d2g1(j,:) = [ig1 ig2 ip1 ip2 d2g1_full(ind_d2g1(j))];
        end
        clear d2g1_full;
    end
end

%% clear variables that are not used any more
clear dys dyy0 dyy0ex0 d2ys
clear rp_static rpp_static
clear g1_static g1p_static g1p g1pp
clear g2_static g2p tmp_g2 g3_tmp
clear ind_d2g1 ind_d2g1_tmp ind_part4 i j i1 i2 ig1 ig2 I II JJ ip1 ip2 is
if order == 1
    clear g2 g3
elseif order == 2
    clear g3
end

%% Construct Jacobian (wrt all selected parameters) of Sigma_e and Corr_e for Gaussian innovations
dSigma_e = zeros(exo_nbr,exo_nbr,totparam_nbr); %initialize
dCorrelation_matrix = zeros(exo_nbr,exo_nbr,totparam_nbr);  %initialize
% Compute Jacobians wrt stderr parameters (these come first)
% note that derivatives wrt stderr parameters are zero by construction for Corr_e
if ~isempty(indpstderr)
    for jp = 1:stderrparam_nbr
        dSigma_e(indpstderr(jp),indpstderr(jp),jp) = 2*sqrt(Sigma_e(indpstderr(jp),indpstderr(jp)));
        if isdiag(Sigma_e) == 0 % if there are correlated errors add cross derivatives
            indotherex0 = 1:exo_nbr;
            indotherex0(indpstderr(jp)) = [];
            for kk = indotherex0
                dSigma_e(indpstderr(jp), kk, jp) = Correlation_matrix(indpstderr(jp),kk)*sqrt(Sigma_e(kk,kk));
                dSigma_e(kk, indpstderr(jp), jp) = dSigma_e(indpstderr(jp), kk, jp); %symmetry
            end
        end
    end
end
% Compute Jacobians wrt corr parameters (these come second)
if ~isempty(indpcorr)
    for jp = 1:corrparam_nbr
        dSigma_e(indpcorr(jp,1),indpcorr(jp,2),stderrparam_nbr+jp) = sqrt(Sigma_e(indpcorr(jp,1),indpcorr(jp,1)))*sqrt(Sigma_e(indpcorr(jp,2),indpcorr(jp,2)));
        dSigma_e(indpcorr(jp,2),indpcorr(jp,1),stderrparam_nbr+jp) = dSigma_e(indpcorr(jp,1),indpcorr(jp,2),stderrparam_nbr+jp); %symmetry
        dCorrelation_matrix(indpcorr(jp,1),indpcorr(jp,2),stderrparam_nbr+jp) = 1;
        dCorrelation_matrix(indpcorr(jp,2),indpcorr(jp,1),stderrparam_nbr+jp) = 1;%symmetry
    end
end
% note that derivatives wrt model parameters (these come third) are zero by construction for Sigma_e and Corr_e

%% Analytical computation of Jacobian (wrt selected model parameters) of first-order solution matrices ghx, ghu, and of Om=ghu*Sigma_e*ghu'
%Notation:
%  fy_ = g1(:,nonzeros(klag))  in DR order
%  fy0 = g1(:,nonzeros(kcurr)) in DR order
%  fyp = g1(:,nonzeros(klead)) in DR order
if analytic_derivation_mode == 1
    % The following derivations are based on Iskrev (2010) and its online appendix A.
    % Basic idea is to make use of the implicit function theorem.
    % Define Kalman transition matrix KalmanA = [0 ghx 0], where the first zero spans nstatic columns, and the second zero nfwrd columns
    % At first order we have: F = GAM0*KalmanA - GAM1*KalmanA*KalmanA - GAM2 = 0, where GAM0=fy0, GAM1=-fyp, GAM2 = -fy_
    % Note that F is a function of parameters p and KalmanA, which is also a function of p,therefore, F = F(p,KalmanA(p))=0, and hence,
    % dF = Fp + dF_dKalmanA*dKalmanA = 0 or dKalmanA = - Fp/dF_dKalmanA

    % Some auxiliary matrices
    I_endo = speye(endo_nbr);
    KalmanA = [zeros(endo_nbr,nstatic) ghx zeros(endo_nbr,nfwrd)];

    % Reshape to write first dynamic derivatives in the Magnus and Neudecker style, i.e. dvec(X)/dp
    GAM0  = zeros(endo_nbr,endo_nbr);
    GAM0(:,kcurr~=0,:) = g1(:,nonzeros(kcurr));
    dGAM0 = zeros(endo_nbr,endo_nbr,modparam_nbr);
    dGAM0(:,kcurr~=0,:) = dg1(:,nonzeros(kcurr),:);
    dGAM0 = reshape(dGAM0, endo_nbr*endo_nbr, modparam_nbr);

    GAM1  = zeros(endo_nbr,endo_nbr);
    GAM1(:,klead~=0,:) = -g1(:,nonzeros(klead));
    dGAM1 = zeros(endo_nbr,endo_nbr,modparam_nbr);
    dGAM1(:,klead~=0,:) = -dg1(:,nonzeros(klead),:);
    dGAM1 = reshape(dGAM1, endo_nbr*endo_nbr, modparam_nbr);

    dGAM2 = zeros(endo_nbr,endo_nbr,modparam_nbr);
    dGAM2(:,klag~=0,:)  = -dg1(:,nonzeros(klag),:);
    dGAM2 = reshape(dGAM2, endo_nbr*endo_nbr,  modparam_nbr);

    GAM3  = -g1(:,yy0_nbr+1:end);
    dGAM3 = reshape(-dg1(:,yy0_nbr+1:end,:), endo_nbr*exo_nbr, modparam_nbr);

    % Compute dKalmanA via implicit function
    dF_dKalmanAghx = kron(I_endo,GAM0) - kron(KalmanA',GAM1) - kron(I_endo,GAM1*KalmanA); %equation 31 in Appendix A of Iskrev (2010)
    Fp = kron(KalmanA',I_endo)*dGAM0 - kron( (KalmanA')^2,I_endo)*dGAM1 - dGAM2; %equation 32 in Appendix A of Iskrev (2010)
    dKalmanA = -dF_dKalmanAghx\Fp;

    % Compute dBB from expressions 33 in Iskrev (2010) Appendix A
    MM = GAM0-GAM1*KalmanA; %this corresponds to matrix M in Ratto and Iskrev (2012, page 6)
    invMM = MM\eye(endo_nbr);
    dghu = - kron( (invMM*GAM3)' , invMM ) * ( dGAM0 - kron( KalmanA' , I_endo ) * dGAM1 - kron( I_endo , GAM1 ) * dKalmanA ) + kron( speye(exo_nbr), invMM ) * dGAM3;

    % Add derivatives for stderr and corr parameters, which are zero by construction
    dKalmanA = [zeros(endo_nbr^2, stderrparam_nbr+corrparam_nbr) dKalmanA];
    dghu = [zeros(endo_nbr*exo_nbr, stderrparam_nbr+corrparam_nbr) dghu];

    % Compute dOm = dvec(ghu*Sigma_e*ghu') from expressions 34 in Iskrev (2010) Appendix A
    dOm = kron(I_endo,ghu*Sigma_e)*(pruned_SS.commutation(endo_nbr, exo_nbr)*dghu)...
          + kron(ghu,ghu)*reshape(dSigma_e, exo_nbr^2, totparam_nbr) + kron(ghu*Sigma_e,I_endo)*dghu;

    % Put into tensor notation
    dKalmanA = reshape(dKalmanA, endo_nbr, endo_nbr, totparam_nbr);
    dghx = dKalmanA(:, nstatic+(1:nspred), stderrparam_nbr+corrparam_nbr+1:end); %get rid of zeros and focus on modparams only
    dghu = reshape(dghu, endo_nbr, exo_nbr, totparam_nbr);
    dghu = dghu(:,:,stderrparam_nbr+corrparam_nbr+1:end); %focus only on modparams
    dOm  = reshape(dOm, endo_nbr, endo_nbr, totparam_nbr);
    clear dF_dKalmanAghx Fp dGAM0 dGAM1 dGAM2 dGAM3 MM invMM I_endo

elseif (analytic_derivation_mode == 0 || analytic_derivation_mode == -2)
    % Here we make use of more efficient generalized sylvester equations
    % Notation: ghx_ = ghx(idx_states,:), ghx0 = ghx(kcurr~=0,:), ghxp = ghx(klead~=0,:)
    % Note that at first-order we have the following expressions, which are (numerically) zero:
    %   * for ghx: g1*zx = fyp*ghxp*ghx_ + fy0*ghx0 + fy_ = A*ghx + fy_ = 0
    %              Taking the differential yields a generalized sylvester equation to get dghx: A*dghx + B*dghx*ghx_ = -dg1*zx
    %   * for ghu: g1*zu = A*ghu + fu = 0
    %              Taking the differential yields an invertible equation to get dghu: A*dghu = -(dfu + dA*ghu)
    % INITIALIZATIONS
    % Note that A and B are the perturbation matrices (NOT the Kalman transition matrices)!
    A = zeros(endo_nbr,endo_nbr);
    A(:,kcurr~=0) = g1(:,nonzeros(kcurr));
    A(:,idx_states) = A(:,idx_states) + g1(:,nonzeros(klead))*ghx(klead~=0,:);
    B = zeros(endo_nbr,endo_nbr);
    B(:,nstatic+npred+1:end) = g1(:,nonzeros(klead));
    zx = [eye(nspred);
          ghx(kcurr~=0,:);
          ghx(klead~=0,:)*ghx(idx_states,:);
          zeros(exo_nbr,nspred)];
    dRHSx = zeros(endo_nbr,nspred,modparam_nbr);
    for jp=1:modparam_nbr
        dRHSx(:,:,jp) = -dg1(:,kyy0,jp)*zx(1:yy0_nbr,:);
    end
    %use iterated generalized sylvester equation to compute dghx
    dghx = sylvester3(A,B,ghx(idx_states,:),dRHSx);
    flag = 1; icount = 0;
    while flag && icount < 4
        [dghx, flag] = sylvester3a(dghx,A,B,ghx(idx_states,:),dRHSx);
        icount = icount+1;
    end

    %Compute dOm, dghu, dA, dB
    dOm = zeros(endo_nbr,endo_nbr,totparam_nbr); %as Om=ghu*Sigma_e*ghu', we need to use totparam_nbr, because there is also a contribution from stderr and corr parameters, which we compute after modparams
    dghu = zeros(endo_nbr,exo_nbr,modparam_nbr);
    dA = zeros(endo_nbr,endo_nbr,modparam_nbr); %dA is also needed at higher orders
    dA(:,kcurr~=0,:) = dg1(:,nonzeros(kcurr),:);
    invA = inv(A); %also needed at higher orders
    for jp=1:modparam_nbr
        dA(:,idx_states,jp) = dA(:,idx_states,jp) + dg1(:,nonzeros(klead),jp)*ghx(klead~=0,:) + g1(:,nonzeros(klead))*dghx(klead~=0,:,jp);
        dghu(:,:,jp) = -invA*( dg1(:,yy0_nbr+1:end,jp) + dA(:,:,jp)*ghu);
        dOm(:,:,stderrparam_nbr+corrparam_nbr+jp) = dghu(:,:,jp)*Sigma_e*ghu' + ghu*Sigma_e*dghu(:,:,jp)';
    end
    %add stderr and corr derivatives to Om=ghu*Sigma_e*ghu'
    if ~isempty(indpstderr)
        for jp = 1:stderrparam_nbr
            dOm(:,:,jp) = ghu*dSigma_e(:,:,jp)*ghu';
        end
    end
    if ~isempty(indpcorr)
        for jp = 1:corrparam_nbr
            dOm(:,:,stderrparam_nbr+jp) = ghu*dSigma_e(:,:,stderrparam_nbr+jp)*ghu';
        end
    end
end

%% Analytical computation of Jacobian (wrt selected model parameters) of second-order solution matrices ghxx,ghxu,ghuu and ghs2
if order > 1
    % Notation: ghxx_ = ghxx(idx_states,:), ghxx0 = ghxx(kcurr~=0,:), ghxxp = ghxx(klead~=0,:)
    %           and similar for ghxu, ghuu and ghs2
    % Note that at second-order we have the following expressions, which are (numerically) zero:
    %   * for ghxx: A*ghxx + B*ghxx*kron(ghx_,ghx_) + g2*kron(zx,zx) = 0
    %               Taking the differential yields a generalized sylvester equation to get dghxx: A*dghxx + B*dghxx*kron(ghx_,ghx_) = RHSxx
    %   * for ghxu: A*ghxu + B*ghxx*kron(ghx_,ghu_) + g2*kron(zx,zu) = 0
    %               Taking the differential yields an invertible equation to get dghxu: A*dghxu = RHSxu
    %   * for ghuu: A*ghuu + B*ghxx*kron(ghu_,ghu_) + g2*kron(zu,zu) = 0
    %               Taking the differential yields an invertible equation to get dghuu: A*dghuu = RHSuu
    %   * for ghs2: Ahs2*ghs2 + (gg2_{++}*kron(ghup,ghup) + fyp*ghuup)*vec(Sigma_e) = 0
    %               Taking the differential yields an invertible equation to get dghs2: S*dghs2 = RHSs2
    %   * due to certainty equivalence and zero mean shocks, we note that ghxs and ghus are zero, and thus not computed
    dB = zeros(endo_nbr,endo_nbr,modparam_nbr); %this matrix is also needed at higher orders
    dB(:,nstatic+npred+1:end,:) = dg1(:,nonzeros(klead),:);
    S = A + B; %needed for dghs2, but also at higher orders
    dS = dA + dB;
    invS = inv(S);
    G_x_x = kron(ghx(idx_states,:),ghx(idx_states,:));
    dG_x_x = zeros(size(G_x_x,1),size(G_x_x,2),modparam_nbr);
    dzx = zeros(size(zx,1),size(zx,2),modparam_nbr);
    dRHSghxx = zeros(endo_nbr,nspred^2,modparam_nbr);
    for jp=1:modparam_nbr
        dzx(:,:,jp) = [zeros(nspred,nspred);
                       dghx(kcurr~=0,:,jp);
                       dghx(klead~=0,:,jp)*ghx(idx_states,:) + ghx(klead~=0,:)*dghx(idx_states,:,jp);
                       zeros(exo_nbr,nspred)];
        dRHS_1 = sparse_hessian_times_B_kronecker_C(dg2(:,k2yy0ex0(kyy0,kyy0)+(jp-1)*yy0ex0_nbr^2),zx(1:yy0_nbr,:),threads_BC);
        dRHS_2 = sparse_hessian_times_B_kronecker_C(g2(:,k2yy0ex0(kyy0,kyy0)),dzx(1:yy0_nbr,:,jp),zx(1:yy0_nbr,:),threads_BC);
        dRHS_3 = sparse_hessian_times_B_kronecker_C(g2(:,k2yy0ex0(kyy0,kyy0)),zx(1:yy0_nbr,:),dzx(1:yy0_nbr,:,jp),threads_BC);
        dG_x_x(:,:,jp) = kron(dghx(idx_states,:,jp),ghx(idx_states,:)) + kron(ghx(idx_states,:),dghx(idx_states,:,jp));
        dRHSghxx(:,:,jp) = -( (dRHS_1+dRHS_2+dRHS_3) + dA(:,:,jp)*ghxx + dB(:,:,jp)*ghxx*G_x_x + B*ghxx*dG_x_x(:,:,jp) );
    end
    %use iterated generalized sylvester equation to compute dghxx
    dghxx = sylvester3(A,B,G_x_x,dRHSghxx);
    flag = 1; icount = 0;
    while flag && icount < 4
        [dghxx, flag] = sylvester3a(dghxx,A,B,G_x_x,dRHSghxx);
        icount = icount+1;
    end
    zu = [zeros(nspred,exo_nbr);
          ghu(kcurr~=0,:);
          ghx(klead~=0,:)*ghu(idx_states,:);
          eye(exo_nbr)];
    abcOutxu = A_times_B_kronecker_C(ghxx,ghx(idx_states,:),ghu(idx_states,:)); %auxiliary expressions for dghxu
    abcOutuu = A_times_B_kronecker_C(ghxx,ghu(idx_states,:));                   %auxiliary expressions for dghuu
    RHSs2 = sparse_hessian_times_B_kronecker_C(g2(:,k2yy0ex0(nonzeros(klead),nonzeros(klead))), ghu(klead~=0,:),threads_BC);
    RHSs2 = RHSs2 + g1(:,nonzeros(klead))*ghuu(klead~=0,:);
    dzu = zeros(size(zu,1),size(zu,2),modparam_nbr);
    dghxu = zeros(endo_nbr,nspred*exo_nbr,modparam_nbr);
    dghuu = zeros(endo_nbr,exo_nbr*exo_nbr,modparam_nbr);
    dghs2 = zeros(endo_nbr,totparam_nbr); %note that for modparam we ignore the contribution by dSigma_e and add it after the computations for stderr and corr
    for jp=1:modparam_nbr
        dzu(:,:,jp) = [zeros(nspred,exo_nbr);
                       dghu(kcurr~=0,:,jp);
                       dghx(klead~=0,:,jp)*ghu(idx_states,:) + ghx(klead~=0,:)*dghu(idx_states,:,jp);
                       zeros(exo_nbr,exo_nbr)];
        %compute dghxu
        dRHS_1 = sparse_hessian_times_B_kronecker_C(dg2(:,k2yy0ex0(kyy0ex0,kyy0ex0)+(jp-1)*yy0ex0_nbr^2),zx,zu,threads_BC);
        dRHS_2 = sparse_hessian_times_B_kronecker_C(g2(:,k2yy0ex0(kyy0ex0,kyy0ex0)),dzx(:,:,jp),zu,threads_BC);
        dRHS_3 = sparse_hessian_times_B_kronecker_C(g2(:,k2yy0ex0(kyy0ex0,kyy0ex0)),zx,dzu(:,:,jp),threads_BC);
        dabcOut_1 = A_times_B_kronecker_C(dghxx(:,:,jp),ghx(idx_states,:),ghu(idx_states,:));
        dabcOut_2 = A_times_B_kronecker_C(ghxx,dghx(idx_states,:,jp),ghu(idx_states,:));
        dabcOut_3 = A_times_B_kronecker_C(ghxx,ghx(idx_states,:),dghu(idx_states,:,jp));
        dghxu(:,:,jp) = -invA * ( dA(:,:,jp)*ghxu + (dRHS_1+dRHS_2+dRHS_3) + dB(:,:,jp)*abcOutxu + B*(dabcOut_1+dabcOut_2+dabcOut_3) );
        %compute dghuu
        dRHS_1 = sparse_hessian_times_B_kronecker_C(dg2(:,k2yy0ex0(kyy0ex0,kyy0ex0)+(jp-1)*yy0ex0_nbr^2),zu,threads_BC);
        dRHS_2 = sparse_hessian_times_B_kronecker_C(g2(:,k2yy0ex0(kyy0ex0,kyy0ex0)),dzu(:,:,jp),zu,threads_BC);
        dRHS_3 = sparse_hessian_times_B_kronecker_C(g2(:,k2yy0ex0(kyy0ex0,kyy0ex0)),zu,dzu(:,:,jp),threads_BC);
        dabcOut_1 = A_times_B_kronecker_C(dghxx(:,:,jp),ghu(idx_states,:));
        dabcOut_2 = A_times_B_kronecker_C(ghxx,dghu(idx_states,:,jp),ghu(idx_states,:));
        dabcOut_3 = A_times_B_kronecker_C(ghxx,ghu(idx_states,:),dghu(idx_states,:,jp));
        dghuu(:,:,jp) = -invA * ( dA(:,:,jp)*ghuu + (dRHS_1+dRHS_2+dRHS_3) + dB(:,:,jp)*abcOutuu + B*(dabcOut_1+dabcOut_2+dabcOut_3) );
        %compute dghs2
        dRHS_1 = sparse_hessian_times_B_kronecker_C(dg2(:,k2yy0ex0(nonzeros(klead),nonzeros(klead))+(jp-1)*yy0ex0_nbr^2), ghu(klead~=0,:),threads_BC);
        dRHS_2 = sparse_hessian_times_B_kronecker_C(g2(:,k2yy0ex0(nonzeros(klead),nonzeros(klead))), dghu(klead~=0,:,jp), ghu(klead~=0,:),threads_BC);
        dRHS_3 = sparse_hessian_times_B_kronecker_C(g2(:,k2yy0ex0(nonzeros(klead),nonzeros(klead))), ghu(klead~=0,:), dghu(klead~=0,:,jp),threads_BC);
        dghs2(:,stderrparam_nbr+corrparam_nbr+jp) = -invS * ( dS(:,:,jp)*ghs2 + ( (dRHS_1+dRHS_2+dRHS_3) + dg1(:,nonzeros(klead),jp)*ghuu(klead~=0,:) + g1(:,nonzeros(klead))*dghuu(klead~=0,:,jp) )*Sigma_e(:) );
    end
    %add contributions by dSigma_e to dghs2
    if ~isempty(indpstderr)
        for jp = 1:stderrparam_nbr
            dghs2(:,jp) = -invS * RHSs2*vec(dSigma_e(:,:,jp));
        end
    end
    if ~isempty(indpcorr)
        for jp = 1:corrparam_nbr
            dghs2(:,stderrparam_nbr+jp) = -invS * RHSs2*vec(dSigma_e(:,:,stderrparam_nbr+jp));
        end
    end
end

if order > 2
    % NOTE: The computations can be improved significantly, particularly for larger models [to do: @wmutschl]
    % Notation: ghxxx_ = ghxxx(idx_states,:), ghxxx0 = ghxxx(kcurr~=0,:), ghxxxp = ghxxx(klead~=0,:)
    %           and similar for ghxxu, ghxuu, ghuuu, ghxss, ghuss

    % Note that at third-order we have the following expressions, which are (numerically) zero, given suitable tensor-unfolding permuation matrices P:
    %   * for ghxxx: A*ghxxx + B*ghxxx*kron(ghx_,kron(ghx_,ghx_)) + g3*kron(kron(zx,zx),zx) + g2*kron(zx,zxx)*P_x_xx + fyp*ghxxp*kron(ghx_,ghxx_)*P_x_xx = 0
    %                Taking the differential yields a generalized sylvester equation to get dghxxx: A*dghxxx + B*dghxxx*kron(ghx_,kron(ghx_,ghx_)) = RHSxxx
    %   * for ghxxu: A*ghxxu + B*ghxxx*kron(ghx_,kron(ghx_,ghu_)) + gg3*kron(zx,kron(zx,zu)) + gg2*(kron(zx,zxu)*P_x_xu + kron(zxx,zu)) + fyp*ghxxp*(kron(ghx_,ghxu_)*P_x_xu + kron(ghxx_,ghu_)) = 0
    %                Taking the differential yields an invertible equation to get dghxxu: A*dghxxu = RHSxxu
    %   * for ghxuu: A*ghxuu + B*ghxxx*kron(ghx_,kron(ghu_,ghu_)) + gg3*kron(zx,kron(zu,zu)) + gg2*(kron(zxu,zu)*P_xu_u + kron(zx,zuu)) + fyp*ghxxp*(kron(ghxu_,ghu_)*Pxu_u + kron(ghx_,ghuu_)) = 0
    %                Taking the differential yields an invertible equation to get dghxuu: A*dghxuu = RHSxuu
    %   * for ghuuu: A*ghuuu + B*ghxxx*kron(ghu_,kron(ghu_,ghu_)) + gg3*kron(kron(zu,zu),zu) + gg2*kron(zu,zuu)*P_u_uu + fyp*ghxxp*kron(ghu_,ghuu_)*P_u_uu = 0
    %                Taking the differential yields an invertible equation to get dghuuu: A*dghuuu = RHSuuu
    %   * for ghxss: A*ghxss + B*ghxss*ghx_ + fyp*ghxxp*kron(ghx_,ghss_) + gg2*kron(zx,zss) + Fxupup*kron(Ix,Sigma_e(:)) = 0
    %                Taking the differential yields a generalized sylvester equation to get dghxss: A*dghxss + B*dghxss*ghx_ = RHSxss
    %   * for ghuss: A*ghuss + B*ghxss*ghu_ + gg2*kron(zu,zss) + fyp*ghxxp*kron(ghu_,ghss_) + Fuupup*kron(Iu,Sigma_e(:)) = 0
    %                Taking the differential yields an invertible equation to get dghuss: A*dghuss = RHSuss
    %   * due to certainty equivalence and Gaussian shocks, we note that ghxxs, ghxus, ghuus, and ghsss are zero and thus not computed

    % permutation matrices
    id_xxx = reshape(1:nspred^3,1,nspred,nspred,nspred);
    id_uux = reshape(1:nspred*exo_nbr^2,1,exo_nbr,exo_nbr,nspred);
    id_uxx = reshape(1:nspred^2*exo_nbr,1,exo_nbr,nspred,nspred);
    id_uuu = reshape(1:exo_nbr^3,1,exo_nbr,exo_nbr,exo_nbr);
    I_xxx = speye(nspred^3);
    I_xxu = speye(nspred^2*exo_nbr);
    I_xuu = speye(nspred*exo_nbr^2);
    I_uuu = speye(exo_nbr^3);
    P_x_xx = I_xxx(:,ipermute(id_xxx,[1,3,4,2])) + I_xxx(:,ipermute(id_xxx,[1,2,4,3])) + I_xxx(:,ipermute(id_xxx,[1,2,3,4]));
    P_x_xu = I_xxu(:,ipermute(id_uxx,[1,2,3,4])) + I_xxu(:,ipermute(id_uxx,[1,2,4,3]));
    P_xu_u = I_xuu(:,ipermute(id_uux,[1,2,3,4])) + I_xuu(:,ipermute(id_uux,[1,3,2,4]));
    P_u_uu = I_uuu(:,ipermute(id_uuu,[1,3,4,2])) + I_uuu(:,ipermute(id_uuu,[1,2,4,3])) + I_uuu(:,ipermute(id_uuu,[1,2,3,4]));
    P_uu_u = I_uuu(:,ipermute(id_uuu,[1,2,3,4])) + I_uuu(:,ipermute(id_uuu,[1,3,4,2]));

    zxx = [spalloc(nspred,nspred^2,0);
           ghxx(kcurr~=0,:);
           ghxx(klead~=0,:)*G_x_x + ghx(klead~=0,:)*ghxx(idx_states,:);
           spalloc(exo_nbr,nspred^2,0)];
    G_x_x_x = kron(ghx(idx_states,:), G_x_x);
    G_x_xx = kron(ghx(idx_states,:),ghxx(idx_states,:));
    Z_x_x = kron(zx,zx);
    Z_x_x_x = kron(zx,Z_x_x);
    Z_x_xx = kron(zx,zxx);
    fyp_ghxxp = sparse(g1(:,nonzeros(klead))*ghxx(klead~=0,:));
    B_ghxxx = B*ghxxx;
    dzxx = zeros(size(zxx,1),size(zxx,2),modparam_nbr);
    dfyp_ghxxp = zeros(size(fyp_ghxxp,1),size(fyp_ghxxp,2),modparam_nbr);

    dRHSghxxx = zeros(endo_nbr,nspred^3,modparam_nbr);
    for jp=1:modparam_nbr
        dzxx(:,:,jp) = [zeros(nspred,nspred^2);
                        dghxx(kcurr~=0,:,jp);
                        dghxx(klead~=0,:,jp)*G_x_x + ghxx(klead~=0,:)*dG_x_x(:,:,jp) + dghx(klead~=0,:,jp)*ghxx(idx_states,:) + ghx(klead~=0,:)*dghxx(idx_states,:,jp);
                        zeros(exo_nbr,nspred^2)];
        dG_x_x_x = kron(dghx(idx_states,:,jp),G_x_x) + kron(ghx(idx_states,:),dG_x_x(:,:,jp));
        dG_x_xx = kron(dghx(idx_states,:,jp),ghxx(idx_states,:)) + kron(ghx(idx_states,:),dghxx(idx_states,:,jp));
        dZ_x_x = kron(dzx(:,:,jp), zx) + kron(zx, dzx(:,:,jp));
        dZ_x_x_x = kron(dzx(:,:,jp), Z_x_x) + kron(zx, dZ_x_x);
        dZ_x_xx = kron(dzx(:,:,jp), zxx) + kron(zx, dzxx(:,:,jp));
        dfyp_ghxxp(:,:,jp) = dg1(:,nonzeros(klead),jp)*ghxx(klead~=0,:) + g1(:,nonzeros(klead))*dghxx(klead~=0,:,jp);
        dRHSghxxx(:,:,jp) = dA(:,:,jp)*ghxxx + dB(:,:,jp)*ghxxx*G_x_x_x + B_ghxxx*dG_x_x_x;
        dRHSghxxx(:,:,jp) = dRHSghxxx(:,:,jp) + dg3(:,k3yy0ex0(kyy0ex0,kyy0ex0,kyy0ex0)+(jp-1)*yy0ex0_nbr^3)*Z_x_x_x + g3(:,k3yy0ex0(kyy0ex0,kyy0ex0,kyy0ex0))*dZ_x_x_x;
        dRHSghxxx(:,:,jp) = dRHSghxxx(:,:,jp) + dg2(:,k2yy0ex0(kyy0ex0,kyy0ex0)+(jp-1)*yy0ex0_nbr^2)*Z_x_xx*P_x_xx + g2(:,k2yy0ex0(kyy0ex0,kyy0ex0))*dZ_x_xx*P_x_xx;
        dRHSghxxx(:,:,jp) = dRHSghxxx(:,:,jp) + dfyp_ghxxp(:,:,jp)*G_x_xx*P_x_xx + fyp_ghxxp*dG_x_xx*P_x_xx;
    end
    dRHSghxxx = -dRHSghxxx;
    %use iterated generalized sylvester equation to compute dghxxx
    dghxxx = sylvester3(A,B,G_x_x_x,dRHSghxxx);
    flag = 1; icount = 0;
    while flag && icount < 4
        [dghxxx, flag] = sylvester3a(dghxxx,A,B,G_x_x_x,dRHSghxxx);
        icount = icount+1;
    end

    %Auxiliary expressions for dghxxu, dghxuu, dghuuu, dghxss, dghuss
    G_x_u = kron(ghx(idx_states,:),ghu(idx_states,:));
    G_u_u = kron(ghu(idx_states,:),ghu(idx_states,:));
    zxu = [zeros(nspred,nspred*exo_nbr);
           ghxu(kcurr~=0,:);
           ghxx(klead~=0,:)*G_x_u + ghx(klead~=0,:)*ghxu(idx_states,:);
           zeros(exo_nbr,exo_nbr*nspred)];
    zuu = [zeros(nspred,exo_nbr^2);
           ghuu(kcurr~=0,:);
           ghxx(klead~=0,:)*G_u_u + ghx(klead~=0,:)*ghuu(idx_states,:);
           zeros(exo_nbr,exo_nbr^2)];
    Z_x_u = kron(zx,zu);
    Z_u_u = kron(zu,zu);
    Z_x_xu = kron(zx,zxu);
    Z_xx_u = kron(zxx,zu);
    Z_xu_u = kron(zxu,zu);
    Z_x_uu = kron(zx,zuu);
    Z_u_uu = kron(zu,zuu);
    Z_x_x_u = kron(Z_x_x,zu);
    Z_x_u_u = kron(Z_x_u,zu);
    Z_u_u_u = kron(Z_u_u,zu);
    G_x_xu = kron(ghx(idx_states,:),ghxu(idx_states,:));
    G_xx_u = kron(ghxx(idx_states,:),ghu(idx_states,:));
    G_xu_u = kron(ghxu(idx_states,:),ghu(idx_states,:));
    G_x_uu = kron(ghx(idx_states,:),ghuu(idx_states,:));
    G_u_uu = kron(ghu(idx_states,:),ghuu(idx_states,:));
    G_x_x_u = kron(G_x_x,ghu(idx_states,:));
    G_x_u_u = kron(G_x_u,ghu(idx_states,:));
    G_u_u_u = kron(G_u_u,ghu(idx_states,:));
    aux_ZP_x_xu_Z_xx_u = Z_x_xu*P_x_xu + Z_xx_u;
    aux_ZP_xu_u_Z_x_uu = Z_xu_u*P_xu_u + Z_x_uu;
    aux_GP_x_xu_G_xx_u = G_x_xu*P_x_xu + G_xx_u;
    aux_GP_xu_u_G_x_uu = G_xu_u*P_xu_u + G_x_uu;
    dghxxu = zeros(endo_nbr,nspred^2*exo_nbr,modparam_nbr);
    dghxuu = zeros(endo_nbr,nspred*exo_nbr^2,modparam_nbr);
    dghuuu = zeros(endo_nbr,exo_nbr^3,modparam_nbr);

    %stuff for ghxss
    zup = [zeros(nspred,exo_nbr);
           zeros(length(nonzeros(kcurr)),exo_nbr);
           ghu(klead~=0,:);
           zeros(exo_nbr,exo_nbr)];
    zss = [zeros(nspred,1);
           ghs2(kcurr~=0,:);
           ghs2(klead~=0,:) + ghx(klead~=0,:)*ghs2(idx_states,:);
           zeros(exo_nbr,1)];
    zxup = [zeros(nspred,nspred*exo_nbr);
            zeros(length(nonzeros(kcurr)),nspred*exo_nbr);
            ghxu(klead~=0,:)*kron(ghx(idx_states,:),eye(exo_nbr));
            zeros(exo_nbr,nspred*exo_nbr)];
    zupup = [zeros(nspred,exo_nbr^2);
             zeros(length(nonzeros(kcurr)),exo_nbr^2);
             ghuu(klead~=0,:);
             zeros(exo_nbr,exo_nbr^2)];
    G_x_ss = kron(ghx(idx_states,:),ghs2(idx_states,:));
    Z_x_ss = kron(zx,zss);
    Z_up_up = kron(zup,zup);
    Z_xup_up = kron(zxup,zup);
    Z_x_upup = kron(zx,zupup);
    Z_x_up_up = kron(zx,Z_up_up);
    aux_ZP_xup_up_Z_x_upup = Z_xup_up*P_xu_u + Z_x_upup;
    Fxupup = g3(:,k3yy0ex0(kyy0ex0,kyy0ex0,kyy0ex0))*Z_x_up_up + g2(:,k2yy0ex0(kyy0ex0,kyy0ex0))*aux_ZP_xup_up_Z_x_upup + g1(:,nonzeros(klead))*ghxuu(klead~=0,:)*kron(ghx(idx_states,:),eye(exo_nbr^2));
    Ix_vecSig_e = kron(speye(nspred),Sigma_e(:));
    dRHSxss = zeros(endo_nbr,nspred,totparam_nbr);

    %stuff for ghuss
    zuup = [zeros(nspred,exo_nbr^2);
            zeros(length(nonzeros(kcurr)),exo_nbr^2);
            ghxu(klead~=0,:)*kron(ghu(idx_states,:),eye(exo_nbr));
            zeros(exo_nbr,exo_nbr^2)];
    G_u_ss = kron(ghu(idx_states,:),ghs2(idx_states,:));
    Z_u_ss = kron(zu,zss);
    Z_u_upup = kron(zu,zupup);
    Z_uup_up = kron(zuup,zup);
    Z_u_up_up = kron(zu,Z_up_up);
    aux_ZP_uup_up_Z_u_upup = Z_uup_up*P_uu_u + Z_u_upup;
    Fuupup = g3(:,k3yy0ex0(kyy0ex0,kyy0ex0,kyy0ex0))*Z_u_up_up + g2(:,k2yy0ex0(kyy0ex0,kyy0ex0))*aux_ZP_uup_up_Z_u_upup + g1(:,nonzeros(klead))*ghxuu(klead~=0,:)*kron(ghu(idx_states,:),eye(exo_nbr^2));
    Iu_vecSig_e = kron(speye(exo_nbr),Sigma_e(:));
    dRHSuss = zeros(endo_nbr,exo_nbr,totparam_nbr);

    for jp=1:modparam_nbr
        dG_x_u = kron(dghx(idx_states,:,jp), ghu(idx_states,:)) + kron(ghx(idx_states,:), dghu(idx_states,:,jp));
        dG_u_u = kron(dghu(idx_states,:,jp), ghu(idx_states,:)) + kron(ghu(idx_states,:), dghu(idx_states,:,jp));
        dzxu = [zeros(nspred,nspred*exo_nbr);
                dghxu(kcurr~=0,:,jp);
                dghxx(klead~=0,:,jp)*G_x_u + ghxx(klead~=0,:)*dG_x_u + dghx(klead~=0,:,jp)*ghxu(idx_states,:) + ghx(klead~=0,:)*dghxu(idx_states,:,jp);
                zeros(exo_nbr,nspred*exo_nbr)];
        dzuu = [zeros(nspred,exo_nbr^2);
                dghuu(kcurr~=0,:,jp);
                dghxx(klead~=0,:,jp)*G_u_u + ghxx(klead~=0,:)*dG_u_u + dghx(klead~=0,:,jp)*ghuu(idx_states,:) + ghx(klead~=0,:)*dghuu(idx_states,:,jp);
                zeros(exo_nbr,exo_nbr^2)];
        dG_x_xu = kron(dghx(idx_states,:,jp),ghxu(idx_states,:)) + kron(ghx(idx_states,:),dghxu(idx_states,:,jp));
        dG_x_uu = kron(dghx(idx_states,:,jp),ghuu(idx_states,:)) + kron(ghx(idx_states,:),dghuu(idx_states,:,jp));
        dG_u_uu = kron(dghu(idx_states,:,jp),ghuu(idx_states,:)) + kron(ghu(idx_states,:),dghuu(idx_states,:,jp));
        dG_xx_u = kron(dghxx(idx_states,:,jp),ghu(idx_states,:)) + kron(ghxx(idx_states,:),dghu(idx_states,:,jp));
        dG_xu_u = kron(dghxu(idx_states,:,jp),ghu(idx_states,:)) + kron(ghxu(idx_states,:),dghu(idx_states,:,jp));
        dG_x_x_u = kron(G_x_x,dghu(idx_states,:,jp)) + kron(dG_x_x(:,:,jp),ghu(idx_states,:));
        dG_x_u_u = kron(G_x_u,dghu(idx_states,:,jp)) + kron(dG_x_u,ghu(idx_states,:));
        dG_u_u_u = kron(G_u_u,dghu(idx_states,:,jp)) + kron(dG_u_u,ghu(idx_states,:));
        dZ_x_u = kron(dzx(:,:,jp),zu) + kron(zx,dzu(:,:,jp));
        dZ_u_u = kron(dzu(:,:,jp),zu) + kron(zu,dzu(:,:,jp));
        dZ_x_x_u = kron(dzx(:,:,jp), Z_x_u) + kron(zx, dZ_x_u);
        dZ_x_u_u = kron(dZ_x_u, zu) + kron(Z_x_u, dzu(:,:,jp));
        dZ_u_u_u = kron(dZ_u_u, zu) + kron(Z_u_u, dzu(:,:,jp));
        dZ_xx_u = kron(dzxx(:,:,jp), zu) + kron(zxx, dzu(:,:,jp));
        dZ_xu_u = kron(dzxu, zu) + kron(zxu, dzu(:,:,jp));
        dZ_x_xu = kron(dzx(:,:,jp), zxu) + kron(zx, dzxu);
        dZ_x_uu = kron(dzx(:,:,jp), zuu) + kron(zx, dzuu);
        dZ_u_uu = kron(dzu(:,:,jp), zuu) + kron(zu, dzuu);
        dB_ghxxx = dB(:,:,jp)*ghxxx + B*dghxxx(:,:,jp);
        %Compute dghxxu
        dRHS = dA(:,:,jp)*ghxxu + dB_ghxxx*G_x_x_u + B_ghxxx*dG_x_x_u;
        dRHS = dRHS + dg3(:,k3yy0ex0(kyy0ex0,kyy0ex0,kyy0ex0)+(jp-1)*yy0ex0_nbr^3)*Z_x_x_u + g3(:,k3yy0ex0(kyy0ex0,kyy0ex0,kyy0ex0))*dZ_x_x_u;
        dRHS = dRHS + dg2(:,k2yy0ex0(kyy0ex0,kyy0ex0)+(jp-1)*yy0ex0_nbr^2)*aux_ZP_x_xu_Z_xx_u + g2(:,k2yy0ex0(kyy0ex0,kyy0ex0))*( dZ_x_xu*P_x_xu + dZ_xx_u );
        dRHS = dRHS + dfyp_ghxxp(:,:,jp)*aux_GP_x_xu_G_xx_u + fyp_ghxxp*( dG_x_xu*P_x_xu + dG_xx_u );
        dghxxu(:,:,jp) = invA* (-dRHS);
        %Compute dghxuu
        dRHS = dA(:,:,jp)*ghxuu + dB_ghxxx*G_x_u_u + B_ghxxx*dG_x_u_u;
        dRHS = dRHS + dg3(:,k3yy0ex0(kyy0ex0,kyy0ex0,kyy0ex0)+(jp-1)*yy0ex0_nbr^3)*Z_x_u_u + g3(:,k3yy0ex0(kyy0ex0,kyy0ex0,kyy0ex0))*dZ_x_u_u;
        dRHS = dRHS + dg2(:,k2yy0ex0(kyy0ex0,kyy0ex0)+(jp-1)*yy0ex0_nbr^2)*aux_ZP_xu_u_Z_x_uu + g2(:,k2yy0ex0(kyy0ex0,kyy0ex0))*( dZ_xu_u*P_xu_u + dZ_x_uu );
        dRHS = dRHS + dfyp_ghxxp(:,:,jp)*aux_GP_xu_u_G_x_uu + fyp_ghxxp*( dG_xu_u*P_xu_u + dG_x_uu );
        dghxuu(:,:,jp) = invA* (-dRHS);
        %Compute dghuuu
        dRHS = dA(:,:,jp)*ghuuu + dB_ghxxx*G_u_u_u + B_ghxxx*dG_u_u_u;
        dRHS = dRHS + dg3(:,k3yy0ex0(kyy0ex0,kyy0ex0,kyy0ex0)+(jp-1)*yy0ex0_nbr^3)*Z_u_u_u + g3(:,k3yy0ex0(kyy0ex0,kyy0ex0,kyy0ex0))*dZ_u_u_u;
        dRHS = dRHS + dg2(:,k2yy0ex0(kyy0ex0,kyy0ex0)+(jp-1)*yy0ex0_nbr^2)*Z_u_uu*P_u_uu + g2(:,k2yy0ex0(kyy0ex0,kyy0ex0))*dZ_u_uu*P_u_uu;
        dRHS = dRHS + dfyp_ghxxp(:,:,jp)*G_u_uu*P_u_uu + fyp_ghxxp*dG_u_uu*P_u_uu;
        dghuuu(:,:,jp) = invA* (-dRHS);
        %Compute dRHSxss
        dzup = [zeros(nspred,exo_nbr);
                zeros(length(nonzeros(kcurr)),exo_nbr);
                dghu(klead~=0,:,jp);
                zeros(exo_nbr,exo_nbr)];
        dzss = [zeros(nspred,1);
                dghs2(kcurr~=0,stderrparam_nbr+corrparam_nbr+jp);
                dghs2(klead~=0,stderrparam_nbr+corrparam_nbr+jp) + dghx(klead~=0,:,jp)*ghs2(idx_states,:) + ghx(klead~=0,:)*dghs2(idx_states,stderrparam_nbr+corrparam_nbr+jp);
                zeros(exo_nbr,1)];
        dzxup = [zeros(nspred,nspred*exo_nbr);
                 zeros(length(nonzeros(kcurr)),nspred*exo_nbr);
                 dghxu(klead~=0,:,jp)*kron(ghx(idx_states,:),eye(exo_nbr)) + ghxu(klead~=0,:)*kron(dghx(idx_states,:,jp),eye(exo_nbr));
                 zeros(exo_nbr,nspred*exo_nbr)];
        dzupup = [zeros(nspred,exo_nbr^2);
                  zeros(length(nonzeros(kcurr)),exo_nbr^2);
                  dghuu(klead~=0,:,jp);
                  zeros(exo_nbr,exo_nbr^2)];
        dG_x_ss = kron(dghx(idx_states,:,jp),ghs2(idx_states,:)) + kron(ghx(idx_states,:),dghs2(idx_states,stderrparam_nbr+corrparam_nbr+jp));
        dZ_x_ss = kron(dzx(:,:,jp),zss) + kron(zx,dzss);
        dZ_up_up = kron(dzup,zup) + kron(zup,dzup);
        dZ_xup_up = kron(dzxup,zup) + kron(zxup,dzup);
        dZ_x_upup = kron(dzx(:,:,jp),zupup) + kron(zx,dzupup);
        dZ_x_up_up = kron(dzx(:,:,jp),Z_up_up) + kron(zx,dZ_up_up);
        daux_ZP_xup_up_Z_x_upup = dZ_xup_up*P_xu_u + dZ_x_upup;
        dFxupup = dg3(:,k3yy0ex0(kyy0ex0,kyy0ex0,kyy0ex0)+(jp-1)*yy0ex0_nbr^3)*Z_x_up_up + g3(:,k3yy0ex0(kyy0ex0,kyy0ex0,kyy0ex0))*dZ_x_up_up...
                  + dg2(:,k2yy0ex0(kyy0ex0,kyy0ex0)+(jp-1)*yy0ex0_nbr^2)*aux_ZP_xup_up_Z_x_upup + g2(:,k2yy0ex0(kyy0ex0,kyy0ex0))*daux_ZP_xup_up_Z_x_upup...
                  + dg1(:,nonzeros(klead),jp)*ghxuu(klead~=0,:)*kron(ghx(idx_states,:),eye(exo_nbr^2)) + g1(:,nonzeros(klead))*dghxuu(klead~=0,:,jp)*kron(ghx(idx_states,:),eye(exo_nbr^2)) + g1(:,nonzeros(klead))*ghxuu(klead~=0,:)*kron(dghx(idx_states,:,jp),eye(exo_nbr^2));
        dRHSxss(:,:,stderrparam_nbr+corrparam_nbr+jp) = dA(:,:,jp)*ghxss + dB(:,:,jp)*ghxss*ghx(idx_states,:) + B*ghxss*dghx(idx_states,:,jp);
        dRHSxss(:,:,stderrparam_nbr+corrparam_nbr+jp) = dRHSxss(:,:,stderrparam_nbr+corrparam_nbr+jp) + dfyp_ghxxp(:,:,jp)*G_x_ss + fyp_ghxxp*dG_x_ss;
        dRHSxss(:,:,stderrparam_nbr+corrparam_nbr+jp) = dRHSxss(:,:,stderrparam_nbr+corrparam_nbr+jp) + dg2(:,k2yy0ex0(kyy0ex0,kyy0ex0)+(jp-1)*yy0ex0_nbr^2)*Z_x_ss + g2(:,k2yy0ex0(kyy0ex0,kyy0ex0))*dZ_x_ss;
        dRHSxss(:,:,stderrparam_nbr+corrparam_nbr+jp) = dRHSxss(:,:,stderrparam_nbr+corrparam_nbr+jp) + dFxupup*Ix_vecSig_e; %missing contribution by dSigma_e
        %Compute dRHSuss
        dzuup = [zeros(nspred,exo_nbr^2);
                 zeros(length(nonzeros(kcurr)),exo_nbr^2);
                 dghxu(klead~=0,:,jp)*kron(ghu(idx_states,:),eye(exo_nbr)) + ghxu(klead~=0,:)*kron(dghu(idx_states,:,jp),eye(exo_nbr));
                 zeros(exo_nbr,exo_nbr^2)];
        dG_u_ss = kron(dghu(idx_states,:,jp),ghs2(idx_states,:)) + kron(ghu(idx_states,:),dghs2(idx_states,stderrparam_nbr+corrparam_nbr+jp));
        dZ_u_ss = kron(dzu(:,:,jp),zss) + kron(zu,dzss);
        dZ_u_upup = kron(dzu(:,:,jp),zupup) + kron(zu,dzupup);
        dZ_uup_up = kron(dzuup,zup) + kron(zuup,dzup);
        dZ_u_up_up = kron(dzu(:,:,jp),Z_up_up) + kron(zu,dZ_up_up);
        daux_ZP_uup_up_Z_u_upup = dZ_uup_up*P_uu_u + dZ_u_upup;
        dFuupup = dg3(:,k3yy0ex0(kyy0ex0,kyy0ex0,kyy0ex0)+(jp-1)*yy0ex0_nbr^3)*Z_u_up_up + g3(:,k3yy0ex0(kyy0ex0,kyy0ex0,kyy0ex0))*dZ_u_up_up...
                  + dg2(:,k2yy0ex0(kyy0ex0,kyy0ex0)+(jp-1)*yy0ex0_nbr^2)*aux_ZP_uup_up_Z_u_upup + g2(:,k2yy0ex0(kyy0ex0,kyy0ex0))*daux_ZP_uup_up_Z_u_upup...
                  + dg1(:,nonzeros(klead),jp)*ghxuu(klead~=0,:)*kron(ghu(idx_states,:),eye(exo_nbr^2)) + g1(:,nonzeros(klead))*dghxuu(klead~=0,:,jp)*kron(ghu(idx_states,:),eye(exo_nbr^2)) + g1(:,nonzeros(klead))*ghxuu(klead~=0,:)*kron(dghu(idx_states,:,jp),eye(exo_nbr^2));
        dRHSuss(:,:,stderrparam_nbr+corrparam_nbr+jp) = dA(:,:,jp)*ghuss + dB(:,:,jp)*ghxss*ghu(idx_states,:) + B*ghxss*dghu(idx_states,:,jp); %missing dghxss
        dRHSuss(:,:,stderrparam_nbr+corrparam_nbr+jp) = dRHSuss(:,:,stderrparam_nbr+corrparam_nbr+jp) + dfyp_ghxxp(:,:,jp)*G_u_ss + fyp_ghxxp*dG_u_ss;
        dRHSuss(:,:,stderrparam_nbr+corrparam_nbr+jp) = dRHSuss(:,:,stderrparam_nbr+corrparam_nbr+jp) + dg2(:,k2yy0ex0(kyy0ex0,kyy0ex0)+(jp-1)*yy0ex0_nbr^2)*Z_u_ss + g2(:,k2yy0ex0(kyy0ex0,kyy0ex0))*dZ_u_ss;
        dRHSuss(:,:,stderrparam_nbr+corrparam_nbr+jp) = dRHSuss(:,:,stderrparam_nbr+corrparam_nbr+jp) + dFuupup*Iu_vecSig_e; %contribution by dSigma_e only for stderr and corr params
    end
    %Add contribution for stderr and corr params to dRHSxss and dRHSuss
    if ~isempty(indpstderr)
        for jp = 1:stderrparam_nbr
            dzss = [zeros(nspred,1);
                    dghs2(kcurr~=0,jp);
                    dghs2(klead~=0,jp) + ghx(klead~=0,:)*dghs2(idx_states,jp);
                    zeros(exo_nbr,1)];
            dG_x_ss = kron(ghx(idx_states,:),dghs2(idx_states,jp));
            dZ_x_ss = kron(zx,dzss);
            dRHSxss(:,:,jp) = Fxupup*kron(speye(nspred),vec(dSigma_e(:,:,jp)));
            dRHSxss(:,:,jp) = dRHSxss(:,:,jp) + fyp_ghxxp*dG_x_ss;
            dRHSxss(:,:,jp) = dRHSxss(:,:,jp) + g2(:,k2yy0ex0(kyy0ex0,kyy0ex0))*dZ_x_ss;

            dG_u_ss = kron(ghu(idx_states,:),dghs2(idx_states,jp));
            dZ_u_ss = kron(zu,dzss);
            dRHSuss(:,:,jp) = Fuupup*kron(speye(exo_nbr),vec(dSigma_e(:,:,jp)));
            dRHSuss(:,:,jp) = dRHSuss(:,:,jp) + fyp_ghxxp*dG_u_ss;
            dRHSuss(:,:,jp) = dRHSuss(:,:,jp) + g2(:,k2yy0ex0(kyy0ex0,kyy0ex0))*dZ_u_ss;
        end
    end
    if ~isempty(indpcorr)
        for jp = (stderrparam_nbr+1):(stderrparam_nbr+corrparam_nbr)
            dzss = [zeros(nspred,1);
                    dghs2(kcurr~=0,jp);
                    dghs2(klead~=0,jp) + ghx(klead~=0,:)*dghs2(idx_states,jp);
                    zeros(exo_nbr,1)];
            dG_x_ss = kron(ghx(idx_states,:),dghs2(idx_states,jp));
            dZ_x_ss = kron(zx,dzss);
            dRHSxss(:,:,jp) = Fxupup*kron(speye(nspred),vec(dSigma_e(:,:,jp)));
            dRHSxss(:,:,jp) = dRHSxss(:,:,jp) + fyp_ghxxp*dG_x_ss;
            dRHSxss(:,:,jp) = dRHSxss(:,:,jp) + g2(:,k2yy0ex0(kyy0ex0,kyy0ex0))*dZ_x_ss;

            dG_u_ss = kron(ghu(idx_states,:),dghs2(idx_states,jp));
            dZ_u_ss = kron(zu,dzss);
            dRHSuss(:,:,jp) = Fuupup*kron(speye(exo_nbr),vec(dSigma_e(:,:,jp)));
            dRHSuss(:,:,jp) = dRHSuss(:,:,jp) + fyp_ghxxp*dG_u_ss;
            dRHSuss(:,:,jp) = dRHSuss(:,:,jp) + g2(:,k2yy0ex0(kyy0ex0,kyy0ex0))*dZ_u_ss;
        end
    end
    dRHSxss = -dRHSxss;
    %use iterated generalized sylvester equation to compute dghxss
    dghxss = sylvester3(A,B,ghx(idx_states,:),dRHSxss);
    flag = 1; icount = 0;
    while flag && icount < 4
        [dghxss, flag] = sylvester3a(dghxss,A,B,ghx(idx_states,:),dRHSxss);
        icount = icount+1;
    end
    %Add contribution by dghxss to dRHSuss and compute it
    dghuss = zeros(endo_nbr,exo_nbr,totparam_nbr);
    for jp = 1:totparam_nbr
        dRHS = dRHSuss(:,:,jp) + B*dghxss(:,:,jp)*ghu(idx_states,:);
        dghuss(:,:,jp) = invA* (-dRHS);
    end
end

%% Store into structure
DERIVS.dg1 = dg1;
DERIVS.dSigma_e = dSigma_e;
DERIVS.dCorrelation_matrix = dCorrelation_matrix;
DERIVS.dYss = dYss;
DERIVS.dghx = cat(3,zeros(endo_nbr,nspred,stderrparam_nbr+corrparam_nbr), dghx);
DERIVS.dghu = cat(3,zeros(endo_nbr,exo_nbr,stderrparam_nbr+corrparam_nbr), dghu);
DERIVS.dOm  = dOm;
if order > 1
    DERIVS.dg2 = dg2;
    DERIVS.dghxx = cat(3,zeros(endo_nbr,nspred^2,stderrparam_nbr+corrparam_nbr), dghxx);
    DERIVS.dghxu = cat(3,zeros(endo_nbr,nspred*exo_nbr,stderrparam_nbr+corrparam_nbr), dghxu);
    DERIVS.dghuu = cat(3,zeros(endo_nbr,exo_nbr^2,stderrparam_nbr+corrparam_nbr), dghuu);
    DERIVS.dghs2 = dghs2;
end
if order > 2
    DERIVS.dg3 = dg3;
    DERIVS.dghxxx = cat(3,zeros(endo_nbr,nspred^3,stderrparam_nbr+corrparam_nbr), dghxxx);
    DERIVS.dghxxu = cat(3,zeros(endo_nbr,nspred^2*exo_nbr,stderrparam_nbr+corrparam_nbr), dghxxu);
    DERIVS.dghxuu = cat(3,zeros(endo_nbr,nspred*exo_nbr^2,stderrparam_nbr+corrparam_nbr), dghxuu);
    DERIVS.dghuuu = cat(3,zeros(endo_nbr,exo_nbr^3,stderrparam_nbr+corrparam_nbr), dghuuu);
    DERIVS.dghxss = dghxss;
    DERIVS.dghuss = dghuss;
end

%% Construct Hessian (wrt all selected parameters) of ghx, and Om=ghu*Sigma_e*ghu'
if d2flag
    % Construct Hessian (wrt all selected parameters) of Sigma_e
    % note that we only need to focus on (stderr x stderr), (stderr x corr), (corr x stderr) parameters, because derivatives wrt all other second-cross parameters are zero by construction
    d2Sigma_e = zeros(exo_nbr,exo_nbr,totparam_nbr^2); %initialize full matrix, even though we'll reduce it later to unique upper triangular values
    % Compute Hessian of Sigma_e wrt (stderr x stderr) parameters
    if ~isempty(indp2stderrstderr)
        for jp = 1:stderrparam_nbr
            for ip = 1:jp
                if jp == ip %same stderr parameters
                    d2Sigma_e(indpstderr(jp),indpstderr(jp),indp2stderrstderr(ip,jp)) = 2;
                else %different stderr parameters
                    if isdiag(Sigma_e) == 0 % if there are correlated errors
                        d2Sigma_e(indpstderr(jp),indpstderr(ip),indp2stderrstderr(ip,jp)) = Correlation_matrix(indpstderr(jp),indpstderr(ip));
                        d2Sigma_e(indpstderr(ip),indpstderr(jp),indp2stderrstderr(ip,jp)) = Correlation_matrix(indpstderr(jp),indpstderr(ip)); %symmetry
                    end
                end
            end
        end
    end
    % Compute Hessian of Sigma_e wrt (stderr x corr) parameters
    if ~isempty(indp2stderrcorr)
        for jp = 1:stderrparam_nbr
            for ip = 1:corrparam_nbr
                if indpstderr(jp) == indpcorr(ip,1) %if stderr is equal to first index of corr parameter, then derivative is equal to stderr of second index
                    d2Sigma_e(indpstderr(jp),indpcorr(ip,2),indp2stderrcorr(jp,ip)) = sqrt(Sigma_e(indpcorr(ip,2),indpcorr(ip,2)));
                    d2Sigma_e(indpcorr(ip,2),indpstderr(jp),indp2stderrcorr(jp,ip)) = sqrt(Sigma_e(indpcorr(ip,2),indpcorr(ip,2))); % symmetry
                end
                if indpstderr(jp) == indpcorr(ip,2) %if stderr is equal to second index of corr parameter, then derivative is equal to stderr of first index
                    d2Sigma_e(indpstderr(jp),indpcorr(ip,1),indp2stderrcorr(jp,ip)) = sqrt(Sigma_e(indpcorr(ip,1),indpcorr(ip,1)));
                    d2Sigma_e(indpcorr(ip,1),indpstderr(jp),indp2stderrcorr(jp,ip)) = sqrt(Sigma_e(indpcorr(ip,1),indpcorr(ip,1))); % symmetry
                end
            end
        end
    end
    d2Sigma_e = d2Sigma_e(:,:,indp2tottot2); %focus on upper triangular hessian values only

    % Construct nonzero derivatives wrt to t+1, i.e. GAM1=-f_{y^+} in Villemot (2011)
    GAM1  = zeros(endo_nbr,endo_nbr);
    GAM1(:,klead~=0,:) = -g1(:,nonzeros(klead));
    dGAM1  = zeros(endo_nbr,endo_nbr,modparam_nbr);
    dGAM1(:,klead~=0,:) = -dg1(:,nonzeros(klead),:);
    indind = ismember(d2g1(:,2),nonzeros(klead));
    tmp = d2g1(indind,:);
    tmp(:,end)=-tmp(:,end);
    d2GAM1 = tmp;
    indklead = find(klead~=0);
    for j = 1:size(tmp,1)
        inxinx = (nonzeros(klead)==tmp(j,2));
        d2GAM1(j,2) = indklead(inxinx);
    end

    % Construct nonzero derivatives wrt to t, i.e. GAM0=f_{y^0} in Villemot (2011)
    GAM0  = zeros(endo_nbr,endo_nbr);
    GAM0(:,kcurr~=0,:) = g1(:,nonzeros(kcurr));
    dGAM0  = zeros(endo_nbr,endo_nbr,modparam_nbr);
    dGAM0(:,kcurr~=0,:) = dg1(:,nonzeros(kcurr),:);
    indind = ismember(d2g1(:,2),nonzeros(kcurr));
    tmp = d2g1(indind,:);
    d2GAM0 = tmp;
    indkcurr = find(kcurr~=0);
    for j = 1:size(tmp,1)
        inxinx = (nonzeros(kcurr)==tmp(j,2));
        d2GAM0(j,2) = indkcurr(inxinx);
    end

    % Construct nonzero derivatives wrt to t-1, i.e. GAM2=-f_{y^-} in Villemot (2011)
    % GAM2 = zeros(endo_nbr,endo_nbr);
    % GAM2(:,klag~=0)  = -g1(:,nonzeros(klag));
    % dGAM2 = zeros(endo_nbr,endo_nbr,modparam_nbr);
    % dGAM2(:,klag~=0)  = -dg1(:,nonzeros(klag),:);
    indind = ismember(d2g1(:,2),nonzeros(klag));
    tmp = d2g1(indind,:);
    tmp(:,end) = -tmp(:,end);
    d2GAM2 = tmp;
    indklag = find(klag~=0);
    for j = 1:size(tmp,1)
        inxinx = (nonzeros(klag)==tmp(j,2));
        d2GAM2(j,2) = indklag(inxinx);
    end

    % Construct nonzero derivatives wrt to u_t, i.e. GAM3=-f_{u} in Villemot (2011)
    % GAM3  = -g1(:,yy0_nbr+1:end);
    % dGAM3  = -dg1(:,yy0_nbr+1:end,:);
    cols_ex = yy0_nbr+(1:yy0ex0_nbr);
    indind = ismember(d2g1(:,2),cols_ex);
    tmp = d2g1(indind,:);
    tmp(:,end) = -tmp(:,end);
    d2GAM3 = tmp;
    for j = 1:size(tmp,1)
        inxinx = find(cols_ex==tmp(j,2));
        d2GAM3(j,2) = inxinx;
    end

    clear d2g1 tmp

    % Compute Hessian (wrt selected params) of ghx using generalized sylvester equations, see equations 17 and 18 in Ratto and Iskrev (2012)
    % solves MM*d2KalmanA+N*d2KalmanA*P = QQ where d2KalmanA are second order derivatives (wrt model parameters) of KalmanA
    QQ = zeros(endo_nbr,endo_nbr,floor(sqrt(modparam_nbr2)));
    jcount=0;
    cumjcount=0;
    jinx = [];
    x2x=sparse(endo_nbr*endo_nbr,modparam_nbr2);
    dKalmanA = zeros(endo_nbr,endo_nbr,modparam_nbr);
    dKalmanA(:,idx_states,:) = dghx;
    MM = (GAM0-GAM1*KalmanA);
    invMM = inv(MM);
    for i=1:modparam_nbr
        for j=1:i
            elem1 = (get_2nd_deriv(d2GAM0,endo_nbr,endo_nbr,j,i)-get_2nd_deriv(d2GAM1,endo_nbr,endo_nbr,j,i)*KalmanA);
            elem1 = get_2nd_deriv(d2GAM2,endo_nbr,endo_nbr,j,i)-elem1*KalmanA;
            elemj0 = dGAM0(:,:,j)-dGAM1(:,:,j)*KalmanA;
            elemi0 = dGAM0(:,:,i)-dGAM1(:,:,i)*KalmanA;
            elem2 = -elemj0*dKalmanA(:,:,i)-elemi0*dKalmanA(:,:,j);
            elem2 = elem2 + ( dGAM1(:,:,j)*dKalmanA(:,:,i) + dGAM1(:,:,i)*dKalmanA(:,:,j) )*KalmanA;
            elem2 = elem2 + GAM1*( dKalmanA(:,:,i)*dKalmanA(:,:,j) + dKalmanA(:,:,j)*dKalmanA(:,:,i));
            jcount=jcount+1;
            jinx = [jinx; [j i]];
            QQ(:,:,jcount) = elem1+elem2;
            if jcount==floor(sqrt(modparam_nbr2)) || (j*i)==modparam_nbr^2
                if (j*i)==modparam_nbr^2
                    QQ = QQ(:,:,1:jcount);
                end
                xx2=sylvester3(MM,-GAM1,KalmanA,QQ);
                flag=1;
                icount=0;
                while flag && icount<4
                    [xx2, flag]=sylvester3a(xx2,MM,-GAM1,KalmanA,QQ);
                    icount = icount + 1;
                end
                x2x(:,cumjcount+1:cumjcount+jcount)=reshape(xx2,[endo_nbr*endo_nbr jcount]);
                cumjcount=cumjcount+jcount;
                jcount = 0;
                jinx = [];
            end
        end
    end
    clear xx2;
    jcount = 0;
    icount = 0;
    cumjcount = 0;
    MAX_DIM_MAT = 100000000;
    ncol = max(1,floor(MAX_DIM_MAT/(8*endo_nbr*(endo_nbr+1)/2)));
    ncol = min(ncol, totparam_nbr2);
    d2KalmanA = sparse(endo_nbr*endo_nbr,totparam_nbr2);
    d2Om = sparse(endo_nbr*(endo_nbr+1)/2,totparam_nbr2);
    d2KalmanA_tmp = zeros(endo_nbr*endo_nbr,ncol);
    d2Om_tmp = zeros(endo_nbr*(endo_nbr+1)/2,ncol);
    tmpDir = CheckPath('tmp_derivs',dname);
    offset = stderrparam_nbr+corrparam_nbr;
    %     d2B = zeros(m,n,tot_param_nbr,tot_param_nbr);
    for j=1:totparam_nbr
        for i=1:j
            jcount=jcount+1;
            if j<=offset %stderr and corr parameters
                y = KalmanB*d2Sigma_e(:,:,jcount)*KalmanB';
                d2Om_tmp(:,jcount) = dyn_vech(y);
            else %model parameters
                jind = j-offset;
                iind = i-offset;
                if i<=offset
                    y = dghu(:,:,jind)*dSigma_e(:,:,i)*KalmanB'+KalmanB*dSigma_e(:,:,i)*dghu(:,:,jind)';
                    %                     y(abs(y)<1.e-8)=0;
                    d2Om_tmp(:,jcount) = dyn_vech(y);
                else
                    icount=icount+1;
                    dKalmanAj = reshape(x2x(:,icount),[endo_nbr endo_nbr]);
                    %                     x = get_2nd_deriv(x2x,m,m,iind,jind);%xx2(:,:,jcount);
                    elem1 = (get_2nd_deriv(d2GAM0,endo_nbr,endo_nbr,iind,jind)-get_2nd_deriv(d2GAM1,endo_nbr,endo_nbr,iind,jind)*KalmanA);
                    elem1 = elem1 -( dGAM1(:,:,jind)*dKalmanA(:,:,iind) + dGAM1(:,:,iind)*dKalmanA(:,:,jind) );
                    elemj0 = dGAM0(:,:,jind)-dGAM1(:,:,jind)*KalmanA-GAM1*dKalmanA(:,:,jind);
                    elemi0 = dGAM0(:,:,iind)-dGAM1(:,:,iind)*KalmanA-GAM1*dKalmanA(:,:,iind);
                    elem0 = elemj0*dghu(:,:,iind)+elemi0*dghu(:,:,jind);
                    y = invMM * (get_2nd_deriv(d2GAM3,endo_nbr,exo_nbr,iind,jind)-elem0-(elem1-GAM1*dKalmanAj)*KalmanB);
                    %         d2B(:,:,j+length(indexo),i+length(indexo)) = y;
                    %         d2B(:,:,i+length(indexo),j+length(indexo)) = y;
                    y = y*Sigma_e*KalmanB'+KalmanB*Sigma_e*y'+ ...
                        dghu(:,:,jind)*Sigma_e*dghu(:,:,iind)'+dghu(:,:,iind)*Sigma_e*dghu(:,:,jind)';
                    %                     x(abs(x)<1.e-8)=0;
                    d2KalmanA_tmp(:,jcount) = vec(dKalmanAj);
                    %                     y(abs(y)<1.e-8)=0;
                    d2Om_tmp(:,jcount) = dyn_vech(y);
                end
            end
            if jcount==ncol || i*j==totparam_nbr^2
                d2KalmanA(:,cumjcount+1:cumjcount+jcount) = d2KalmanA_tmp(:,1:jcount);
                %         d2KalmanA(:,:,j+length(indexo),i+length(indexo)) = x;
                %         d2KalmanA(:,:,i+length(indexo),j+length(indexo)) = x;
                d2Om(:,cumjcount+1:cumjcount+jcount) = d2Om_tmp(:,1:jcount);
                %         d2Om(:,:,j+length(indexo),i+length(indexo)) = y;
                %         d2Om(:,:,i+length(indexo),j+length(indexo)) = y;
                save([tmpDir filesep 'd2KalmanA_' int2str(cumjcount+1) '_' int2str(cumjcount+jcount) '.mat'],'d2KalmanA')
                save([tmpDir filesep 'd2Om_' int2str(cumjcount+1) '_'  int2str(cumjcount+jcount) '.mat'],'d2Om')
                cumjcount = cumjcount+jcount;
                jcount=0;
                %         d2KalmanA = sparse(m1*m1,tot_param_nbr*(tot_param_nbr+1)/2);
                %         d2Om = sparse(m1*(m1+1)/2,tot_param_nbr*(tot_param_nbr+1)/2);
                d2KalmanA_tmp = zeros(endo_nbr*endo_nbr,ncol);
                d2Om_tmp = zeros(endo_nbr*(endo_nbr+1)/2,ncol);
            end
        end
    end

    %Store into structure
    DERIVS.d2Yss     = d2Yss;
    DERIVS.d2KalmanA = d2KalmanA;
    DERIVS.d2Om      = d2Om;
end

return

%% AUXILIARY FUNCTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%

function g22 = get_2nd_deriv(gpp,m,n,i,j)
% inputs:
% - gpp: [#second_order_Jacobian_terms by 5] double   Hessian matrix (wrt parameters) of a matrix
%                                                              rows: respective derivative term
%                                                              1st column: equation number of the term appearing
%                                                              2nd column: column number of variable in Jacobian
%                                                              3rd column: number of the first parameter in derivative
%                                                              4th column: number of the second parameter in derivative
%                                                              5th column: value of the Hessian term
% - m:    scalar                                     number of equations
% - n:    scalar                                     number of variables
% - i:    scalar                                     number for which first parameter
% - j:    scalar                                     number for which second parameter

g22=zeros(m,n);
is=find(gpp(:,3)==i);
is=is(find(gpp(is,4)==j));

if ~isempty(is)
    g22(sub2ind([m,n],gpp(is,1),gpp(is,2)))=gpp(is,5)';
end
return

function g22 = get_2nd_deriv_mat(gpp,i,j,npar)
% inputs:
% - gpp: [#second_order_Jacobian_terms by 5] double   Hessian matrix of (wrt parameters) of dynamic Jacobian
%                                                              rows: respective derivative term
%                                                              1st column: equation number of the term appearing
%                                                              2nd column: column number of variable in Jacobian of the dynamic model
%                                                              3rd column: number of the first parameter in derivative
%                                                              4th column: number of the second parameter in derivative
%                                                              5th column: value of the Hessian term
% - i:    scalar                                     number for which model equation
% - j:    scalar                                     number for which variable in Jacobian of dynamic model
% - npar: scalar                                     Number of model parameters, i.e. equals param_nbr
%
% output:
% g22: [npar by npar] Hessian matrix (wrt parameters) of Jacobian of dynamic model for equation i
%                                                    rows: first parameter in Hessian
%                                                    columns: second paramater in Hessian

g22=zeros(npar,npar);
is=find(gpp(:,1)==i);
is=is(find(gpp(is,2)==j));

if ~isempty(is)
    g22(sub2ind([npar,npar],gpp(is,3),gpp(is,4)))=gpp(is,5)';
end
return

function r22 = get_all_resid_2nd_derivs(rpp,m,npar)
% inputs:
% - rpp: [#second_order_residual_terms by 4] double   Hessian matrix (wrt paramters) of model equations
%                                                              rows: respective derivative term
%                                                              1st column: equation number of the term appearing
%                                                              2nd column: number of the first parameter in derivative
%                                                              3rd column: number of the second parameter in derivative
%                                                              4th column: value of the Hessian term
% - m:    scalar                                     Number of residuals (or model equations), i.e. equals endo_nbr
% - npar: scalar                                     Number of model parameters, i.e. equals param_nbr
%
% output:
% r22: [endo_nbr by param_nbr by param_nbr] Hessian matrix of model equations with respect to parameters
%                                                              rows: equations in order of declaration
%                                                              1st columns: first parameter number in derivative
%                                                              2nd columns: second parameter in derivative

r22=zeros(m,npar,npar);

for is=1:size(rpp,1)
    % Keep symmetry in hessian, hence 2 and 3 as well as 3 and 2, i.e. d2f/(dp1 dp2) = d2f/(dp2 dp1)
    r22(rpp(is,1),rpp(is,2),rpp(is,3))=rpp(is,4);
    r22(rpp(is,1),rpp(is,3),rpp(is,2))=rpp(is,4);
end

return

function h2 = get_hess_deriv(hp,i,j,m,npar)
% inputs:
% - hp: [#first_order_Hessian_terms by 5] double   Jacobian matrix (wrt paramters) of dynamic Hessian
%                                                              rows: respective derivative term
%                                                              1st column: equation number of the term appearing
%                                                              2nd column: column number of first variable in Hessian of the dynamic model
%                                                              3rd column: column number of second variable in Hessian of the dynamic model
%                                                              4th column: number of the parameter in derivative
%                                                              5th column: value of the Hessian term
% - i:    scalar                                     number for which model equation
% - j:    scalar                                     number for which first variable in Hessian of dynamic model variable
% - m:    scalar                                     Number of dynamic model variables + exogenous vars, i.e. yy0_nbr + exo_nbr
% - npar: scalar                                     Number of model parameters, i.e. equals param_nbr
%
% output:
% h2: [(yy0_nbr + exo_nbr) by param_nbr] Jacobian matrix (wrt parameters) of dynamic Hessian
%                                                              rows: second dynamic or exogenous variables in Hessian of specific model equation of the dynamic model
%                                                              columns: parameters

h2=zeros(m,npar);
is1=find(hp(:,1)==i);
is=is1(find(hp(is1,2)==j));

if ~isempty(is)
    h2(sub2ind([m,npar],hp(is,3),hp(is,4)))=hp(is,5)';
end

return
