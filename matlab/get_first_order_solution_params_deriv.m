function [dA, dB, dSigma_e, dOm, dYss, dg1, d2A, d2Om, d2Yss] = get_first_order_solution_params_deriv(A, B, estim_params, M, oo, options, kronflag, indpmodel, indpstderr, indpcorr, indvar)
%[dA, dB, dSigma_e, dOm, dYss, dg1, d2A, d2Om, d2Yss] = get_first_order_solution_params_deriv(A, B, estim_params, M, oo, options, kronflag, indpmodel, indpstderr, indpcorr, indvar)
% previously getH.m
% -------------------------------------------------------------------------
% Computes first and second derivatives (with respect to parameters) of
%   (1) reduced-form solution (dA, dB, dSigma_e, dOm, d2A, d2Om)
%   (1) steady-state (dYss, d2Yss)
%   (3) Jacobian (wrt to dynamic variables) of dynamic model (dg1)
% Note that the order in the parameter Jacobians is the following:
% first stderr parameters, second corr parameters, third model parameters
% =========================================================================
% INPUTS
%   A:            [endo_nbr by endo_nbr] Transition matrix from Kalman filter
%                   for all endogenous declared variables, in DR order
%   B:            [endo_nbr by exo_nbr]  Transition matrix from Kalman filter
%                   mapping shocks today to endogenous variables today, in DR order
%   estim_params: [structure] storing the estimation information
%   M:            [structure] storing the model information
%   oo:           [structure] storing the reduced-form solution results
%   options:      [structure] storing the options
%   kronflag:     [scalar] method to compute Jacobians (equal to analytic_derivation_mode in options_ident). Default:0
%                   *  0: efficient sylvester equation method to compute
%                      analytical derivatives as in Ratto & Iskrev (2011)
%                   *  1: kronecker products method to compute analytical
%                      derivatives as in Iskrev (2010)
%                   * -1: numerical two-sided finite difference method to
%                      compute numerical derivatives of all output arguments
%                      using function identification_numerical_objective.m
%                      (previously thet2tau.m)
%                   * -2: numerical two-sided finite difference method to 
%                      compute numerically dYss, dg1, d2Yss and d2g1, the other 
%                      output arguments are computed analytically as in kronflag=0
%   indpmodel:    [modparam_nbr by 1] index of estimated parameters in M_.params; 
%                   corresponds to model parameters (no stderr and no corr) 
%                   in estimated_params block; if estimated_params block is
%                   not available, then all model parameters are selected
%   indpstderr:   [stderrparam_nbr by 1] index of estimated standard errors, 
%                   i.e. for all exogenous variables where "stderr" is given 
%                   in the estimated_params block; if estimated_params block
%                   is not available, then all stderr parameters are selected
%   indpcorr:     [corrparam_nbr by 2] matrix of estimated correlations,
%                   i.e. for all exogenous variables where "corr" is given 
%                   in the estimated_params block; if estimated_params block
%                   is not available, then no corr parameters are selected
%   indvar:       [var_nbr by 1] index of considered (or observed) variables
% -------------------------------------------------------------------------
% OUTPUTS
%   dA:         [var_nbr by var_nbr by totparam_nbr] in DR order
%                   Jacobian (wrt to all parameters) of transition matrix A
%   dB:         [var_nbr by exo_nbr by totparam_nbr] in DR order
%                   Jacobian (wrt to all parameters) of transition matrix B
%   dSigma_e:   [exo_nbr by exo_nbr by totparam_nbr] in declaration order
%                   Jacobian (wrt to all paramters) of M_.Sigma_e
%   dOm:        [var_nbr by var_nbr by totparam_nbr] in DR order
%                   Jacobian (wrt to all paramters) of Om = (B*M_.Sigma_e*B')
%   dYss:       [var_nbr by modparam_nbr] in DR order
%                   Jacobian (wrt model parameters only) of steady state
%   dg1:        [endo_nbr by (dynamicvar_nbr + exo_nbr) by modparam_nbr] in DR order 
%                   Jacobian (wrt to model parameters only) of Jacobian of dynamic model
%   d2A:        [var_nbr*var_nbr by totparam_nbr*(totparam_nbr+1)/2] in DR order
%                   Unique entries of Hessian (wrt all parameters) of transition matrix A
%   d2Om:       [var_nbr*(var_nbr+1)/2 by totparam_nbr*(totparam_nbr+1)/2] in DR order
%                   Unique entries of Hessian (wrt all parameters) of Omega
%   d2Yss:      [var_nbr by modparam_nbr by modparam_nbr] in DR order
%                   Unique entries of Hessian (wrt model parameters only) of steady state
% -------------------------------------------------------------------------
% This function is called by 
%   * dsge_likelihood.m
%   * get_identification_jacobians.m (previously getJJ.m)
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
%   * identification_numerical_objective.m (previously thet2tau.m)
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
fname              = M.fname;
dname              = M.dname;
maximum_exo_lag    = M.maximum_exo_lag;
maximum_exo_lead   = M.maximum_exo_lead;
maximum_endo_lag   = M.maximum_endo_lag;
maximum_endo_lead  = M.maximum_endo_lead;
lead_lag_incidence = M.lead_lag_incidence;
[I,~] = find(lead_lag_incidence');   %I is used to select nonzero columns of the Jacobian of endogenous variables in dynamic model files

ys        = oo.dr.ys;                %steady state of endogenous variables in declaration order
yy0       = oo.dr.ys(I);             %steady state of dynamic (endogenous and auxiliary variables) in DR order
ex0       = oo.exo_steady_state';    %steady state of exogenous variables in declaration order
params0   = M.params;                %values at which to evaluate dynamic, static and param_derivs files
Sigma_e0  = M.Sigma_e;               %covariance matrix of exogenous shocks
Corr_e0   = M.Correlation_matrix;    %correlation matrix of exogenous shocks
stderr_e0 = sqrt(diag(Sigma_e0));    %standard errors of exogenous shocks

param_nbr       = M.param_nbr;        %number of all declared model parameters in mod file
modparam_nbr    = length(indpmodel);  %number of model parameters to be used
stderrparam_nbr = length(indpstderr); %number of stderr parameters to be used
corrparam_nbr   = size(indpcorr,1);   %number of stderr parameters to be used
totparam_nbr    = modparam_nbr + stderrparam_nbr + corrparam_nbr; %total number of parameters to be used

if nargout > 6
    modparam_nbr2 = modparam_nbr*(modparam_nbr+1)/2; %number of unique entries of model parameters only in second-order derivative matrix
    totparam_nbr2 = totparam_nbr*(totparam_nbr+1)/2; %number of unique entries of all parameters in second-order derivative matrix
    %get indices of elements in second derivatives of parameters
    indp2tottot = reshape(1:totparam_nbr^2,totparam_nbr,totparam_nbr);
    indp2stderrstderr = indp2tottot(1:stderrparam_nbr , 1:stderrparam_nbr);
    indp2stderrcorr = indp2tottot(1:stderrparam_nbr , stderrparam_nbr+1:stderrparam_nbr+corrparam_nbr);
    indp2modmod = indp2tottot(stderrparam_nbr+corrparam_nbr+1:stderrparam_nbr+corrparam_nbr+modparam_nbr , stderrparam_nbr+corrparam_nbr+1:stderrparam_nbr+corrparam_nbr+modparam_nbr);
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
end
endo_nbr = size(A,1);      %number of all declared endogenous variables
var_nbr  = length(indvar); %number of considered variables
exo_nbr  = size(B,2);      %number of exogenous shocks in model

if kronflag == -1
% numerical two-sided finite difference method using function identification_numerical_objective.m (previously thet2tau.m) for Jacobian (wrt parameters) of A, B, Sig, Om, Yss, and g1
    para0 = get_all_parameters(estim_params, M);  %get all selected parameters in estimated_params block, stderr and corr come first, then model parameters
    if isempty(para0)
        %if there is no estimated_params block, consider all stderr and all model parameters, but no corr parameters
        para0 = [stderr_e0', params0'];
    end
    %Jacobians (wrt paramters) of steady state, solution matrices A and B, as well as Sigma_e for ALL variables [outputflag=0]
    dYssABSige = fjaco('identification_numerical_objective', para0, 0, estim_params, M, oo, options, indpmodel, indpstderr, indpcorr, indvar);
    M.params = params0;              %make sure values are set back
    M.Sigma_e = Sigma_e0;            %make sure values are set back
    M.Correlation_matrix = Corr_e0 ; %make sure values are set back
    % get Jacobians for Yss, A, and B from dYssABSige
    indYss     = 1:var_nbr;
    indA       = (var_nbr+1):(var_nbr+var_nbr^2);
    indB       = (var_nbr+var_nbr^2+1):(var_nbr+var_nbr^2+var_nbr*exo_nbr);
    indSigma_e = (var_nbr+var_nbr^2+var_nbr*exo_nbr+1):(var_nbr+var_nbr^2+var_nbr*exo_nbr+exo_nbr*(exo_nbr+1)/2);
    dYss = dYssABSige(indYss , stderrparam_nbr+corrparam_nbr+1:end); %in tensor notation only wrt model parameters
    dA = reshape(dYssABSige(indA , :) , [var_nbr var_nbr totparam_nbr]); %in tensor notation
    dB = reshape(dYssABSige(indB , :) , [var_nbr exo_nbr totparam_nbr]); %in tensor notation
    
    dOm = zeros(var_nbr,var_nbr,totparam_nbr); %initialize in tensor notation
    dSigma_e = zeros(exo_nbr,exo_nbr,totparam_nbr);  %initialize in tensor notation
    % get Jacobians of Sigma_e and Om wrt stderr parameters
    if ~isempty(indpstderr)
        for jp=1:stderrparam_nbr
            dSigma_e(:,:,jp) = dyn_unvech(dYssABSige(indSigma_e , jp));
            dOm(:,:,jp) = B*dSigma_e(:,:,jp)*B'; %note that derivatives of B wrt stderr parameters are zero by construction
        end
    end
    % get Jacobians of Sigma_e and Om wrt corr parameters
    if ~isempty(indpcorr)
        for jp=1:corrparam_nbr
            dSigma_e(:,:,stderrparam_nbr+jp) = dyn_unvech(dYssABSige(indSigma_e , stderrparam_nbr+jp));
            dOm(:,:,stderrparam_nbr+jp) = B*dSigma_e(:,:,stderrparam_nbr+jp)*B'; %note that derivatives of B wrt corr parameters are zero by construction
        end
    end
    % get Jacobian of Om wrt model parameters
    if ~isempty(indpmodel)
        for jp=1:modparam_nbr
            dOm(:,:,stderrparam_nbr+corrparam_nbr+jp) = dB(:,:,stderrparam_nbr+corrparam_nbr+jp)*Sigma_e0*B' + B*Sigma_e0*dB(:,:,stderrparam_nbr+corrparam_nbr+jp)'; %note that derivatives of Sigma_e wrt model parameters are zero by construction
        end
    end
    
    %Jacobian (wrt model parameters ONLY) of steady state and of Jacobian of all dynamic model equations [outputflag=-1]    
    dYssg1 = fjaco('identification_numerical_objective', params0(indpmodel), -1, estim_params, M, oo, options, indpmodel, [], [], (1:endo_nbr)');
    M.params = params0;              %make sure values are set back
    M.Sigma_e = Sigma_e0;            %make sure values are set back
    M.Correlation_matrix = Corr_e0 ; %make sure values are set back
    dg1 = reshape(dYssg1(endo_nbr+1:end,:),[endo_nbr, length(yy0)+length(ex0), modparam_nbr]); %get rid of steady state and in tensor notation

    if nargout > 6
        %Hessian (wrt paramters) of steady state, solution matrices A and Om [outputflag=-2]
        % note that hessian_sparse does not take symmetry into account, i.e. compare hessian_sparse.m to hessian.m, but focuses already on unique values, which are duplicated below
        d2YssAOm = hessian_sparse('identification_numerical_objective', para0', options.gstep, -2, estim_params, M, oo, options, indpmodel, indpstderr, indpcorr, indvar);
        M.params = params0;              %make sure values are set back
        M.Sigma_e = Sigma_e0;            %make sure values are set back
        M.Correlation_matrix = Corr_e0 ; %make sure values are set back        
        
        d2A = d2YssAOm(indA , indp2tottot2);    %only unique elements
        d2Om = d2YssAOm(indA(end)+1:end , indp2tottot2);    %only unique elements
        d2Yss = zeros(var_nbr,modparam_nbr,modparam_nbr);   %initialize        
        for j = 1:var_nbr
            d2Yss(j,:,:) = dyn_unvech(full(d2YssAOm(j,indp2modmod2))); %full Hessian for d2Yss, note that here we duplicate unique values for model parameters
        end
        clear d2YssAOm
    end
    
    return %[END OF MAIN FUNCTION]!!!!!
end

if kronflag == -2
% numerical two-sided finite difference method to compute numerically 
% dYss, dg1, d2Yss and d2g1, the rest is computed analytically (kronflag=0) below    
    modpara0 = params0(indpmodel); %focus only on model parameters for dYss, d2Yss and dg1
    [~, g1] = feval([fname,'.dynamic'], yy0, ex0, params0, ys, 1);
    %g1 is [endo_nbr by (dynamicvar_nbr+exo_nbr)] first derivative (wrt all endogenous, exogenous and auxiliary variables) of dynamic model equations, i.e. df/d[yy0;ex0], in DR order
    if nargout > 6
        % computation of d2Yss and d2g1, i.e. second derivative (wrt. parameters) of Jacobian (wrt endogenous and auxilary variables) of dynamic model [outputflag = -1]
        % note that hessian_sparse does not take symmetry into account, i.e. compare hessian_sparse.m to hessian.m, but focuses already on unique values, which are duplicated below
        d2Yssg1 = hessian_sparse('identification_numerical_objective', modpara0, options.gstep, -1, estim_params, M, oo, options, indpmodel, [], [], (1:endo_nbr)');        
        M.params = params0;              %make sure values are set back
        M.Sigma_e = Sigma_e0;            %make sure values are set back
        M.Correlation_matrix = Corr_e0 ; %make sure values are set back
        
        d2Yss = reshape(full(d2Yssg1(1:endo_nbr,:)), [endo_nbr modparam_nbr modparam_nbr]); %put into tensor notation
        for j=1:endo_nbr
            d2Yss(j,:,:) = dyn_unvech(dyn_vech(d2Yss(j,:,:))); %add duplicate values to full hessian
        end
        d2g1_full = d2Yssg1(endo_nbr+1:end,:);
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
        clear d2g1_full;
    end
    %Jacobian (wrt parameters) of steady state and Jacobian of dynamic model equations [outputflag=-1]
    dg1 = fjaco('identification_numerical_objective', modpara0, -1, estim_params, M, oo, options, indpmodel, [], [], (1:endo_nbr)');    
    M.params = params0;              %make sure values are set back
    M.Sigma_e = Sigma_e0;            %make sure values are set back
    M.Correlation_matrix = Corr_e0 ; %make sure values are set back
    dYss = dg1(1:endo_nbr , :);
    dg1 = reshape(dg1(endo_nbr+1 : end , :),[endo_nbr, length(yy0)+length(ex0), modparam_nbr]); %get rid of steady state
elseif (kronflag == 0 || kronflag == 1)
% Analytical method to compute dYss, dg1, d2Yss and d2g1
    [~, g1_static] = feval([fname,'.static'], ys, ex0, params0);
    %g1_static is [endo_nbr by endo_nbr] first-derivative (wrt variables) of static model equations f, i.e. df/dys, in declaration order
    rp_static = feval([fname,'.static_params_derivs'], ys, repmat(ex0, maximum_exo_lag+maximum_exo_lead+1), params0);
    %rp_static is [endo_nbr by param_nbr] first-derivative (wrt parameters) of static model equations f, i.e. df/dparams, in declaration order
    dys = -g1_static\rp_static;
    %use implicit function theorem (equation 5 of Ratto and Iskrev (2011) to compute [endo_nbr by param_nbr] first-derivative (wrt parameters) of steady state analytically, note that dys is in declaration order
    d2ys = zeros(length(ys), param_nbr, param_nbr); %initialize in tensor notation
    if nargout > 6
        [~, ~, g2_static] = feval([fname,'.static'], ys, ex0, params0);
        %g2_static is [endo_nbr by endo_nbr^2] second derivative (wrt variables) of static model equations f, i.e. d(df/dys)/dys, in declaration order
        [~, g1, g2, g3] = feval([fname,'.dynamic'], yy0, ex0, params0, ys, 1);
        %g1 is [endo_nbr by (dynamicvar_nbr+exo_nbr)] first derivative (wrt all endogenous, exogenous and auxiliary variables) of dynamic model equations, i.e. df/d[yy0;ex0], in DR order
        %g2 is [endo_nbr by (dynamicvar_nbr+exo_nbr)^2] second derivative (wrt all endogenous, exogenous and auxiliary variables) of dynamic model equations, i.e. d(df/d[yy0;ex0])/d[yy0;ex0], in DR order
        %g3 is [endo_nbr by (dynamicvar_nbr+exo_nbr)^2] third-derivative (wrt all endogenous, exogenous and auxiliary variables) of dynamic model equations, i.e. d(df/d[yy0;ex0])/d[yy0;ex0], in DR order
        [~, gp_static, rpp_static] = feval([fname,'.static_params_derivs'], ys, ex0, params0);
        %gp_static is [endo_nbr by endo_nbr by param_nbr] first derivative (wrt parameters) of first-derivative (wrt variables) of static model equations f, i.e. (df/dys)/dparams, in declaration order
        %rpp_static are nonzero values and corresponding indices of second derivative (wrt parameters) of static model equations f, i.e. d(df/dparams)/dparams, in declaration order
        rpp_static = get_all_resid_2nd_derivs(rpp_static, length(ys), param_nbr); %make full matrix out of nonzero values and corresponding indices
        %rpp_static is [endo_nbr by param_nbr by param_nbr] second derivative (wrt parameters) of static model equations, i.e. d(df/dparams)/dparams, in declaration order
        if isempty(find(g2_static))
            %auxiliary expression on page 8 of Ratto and Iskrev (2011) is zero, i.e. gam = 0
            for j = 1:param_nbr
                %using the implicit function theorem, equation 15 on page 7 of Ratto and Iskrev (2011)
                d2ys(:,:,j) = -g1_static\rpp_static(:,:,j);
                %d2ys is [endo_nbr by param_nbr by param_nbr] second-derivative (wrt parameters) of steady state, i.e. d(dys/dparams)/dparams, in declaration order
            end
        else
            gam = rpp_static*0; %initialize auxiliary expression on page 8 of Ratto and Iskrev (2011)
            for j = 1:endo_nbr
                tmp_gp_static_dys = (squeeze(gp_static(j,:,:))'*dys);
                gam(j,:,:) = transpose(reshape(g2_static(j,:),[endo_nbr endo_nbr])*dys)*dys + tmp_gp_static_dys + tmp_gp_static_dys';
            end
            for j = 1:param_nbr
                %using the implicit function theorem, equation 15 on page 7 of Ratto and Iskrev (2011)
                d2ys(:,:,j) = -g1_static\(rpp_static(:,:,j)+gam(:,:,j));
                %d2ys is [endo_nbr by param_nbr by param_nbr] second-derivative (wrt parameters) of steady state, i.e. d(dys/dparams)/dparams, in declaration order
            end
            clear gp_static g2_static tmp_gp_static_dys gam
        end
    end
    %handling of steady state for nonstationary variables
    if any(any(isnan(dys)))
        [U,T] = schur(g1_static);
        qz_criterium = options.qz_criterium;
        e1 = abs(ordeig(T)) < qz_criterium-1;
        k = sum(e1);       % Number of non stationary variables.
                           % Number of stationary variables: n = length(e1)-k
        [U,T] = ordschur(U,T,e1);
        T = T(k+1:end,k+1:end);
        %using implicit function theorem, equation 5 of Ratto and Iskrev (2011), in declaration order
        dys = -U(:,k+1:end)*(T\U(:,k+1:end)')*rp_static;
        if nargout > 6
            disp('Computation of d2ys for nonstationary variables is not yet correctly handled if g2_static is nonempty, but continue anyways...')
            for j = 1:param_nbr
                %using implicit function theorem, equation 15 of Ratto and Iskrev (2011), in declaration order
                d2ys(:,:,j) = -U(:,k+1:end)*(T\U(:,k+1:end)')*rpp_static(:,:,j); %THIS IS NOT CORRECT, IF g2_static IS NONEMPTY. WE NEED TO ADD GAM [willi]
            end
        end
    end
    if nargout > 6
        [~, gp, ~, gpp, hp] = feval([fname,'.dynamic_params_derivs'], yy0, ex0, params0, ys, 1, dys, d2ys);
        %gp is [endo_nbr by (dynamicvar_nbr + exo_nbr) by param_nbr] first-derivative (wrt parameters) of first-derivative (wrt all endogenous, auxiliary and exogenous variables) of dynamic model equations, i.e. d(df/dvars)/dparam, in DR order
        %gpp are nonzero values and corresponding indices of second-derivative (wrt parameters) of first-derivative (wrt all endogenous, auxiliary and exogenous variables) of dynamic model equations, i.e. d(d(df/dvars)/dparam)/dparam, in DR order
        %hp are nonzero values and corresponding indices of first-derivative (wrt parameters) of second-derivative (wrt all endogenous, auxiliary and exogenous variables) of dynamic model equations, i.e. d(d(df/dvars)/dvars)/dparam, in DR order
        d2Yss = d2ys(oo.dr.order_var,indpmodel,indpmodel);
        %[endo_nbr by mod_param_nbr by mod_param_nbr], i.e. put into DR order and focus only on model parameters
    else
        [~, gp] = feval([fname,'.dynamic_params_derivs'], yy0, repmat(ex0, [maximum_exo_lag+maximum_exo_lead+1,1]), params0, ys, 1, dys, d2ys);
        %gp is [endo_nbr by (dynamicvar_nbr + exo_nbr) by param_nbr] first-derivative (wrt parameters) of first-derivative (wrt all endogenous, auxiliary and exogenous variables) of dynamic model equations, i.e. d(df/dvars)/dparam, in DR order
        [~, g1, g2 ] = feval([fname,'.dynamic'], yy0, repmat(ex0, [maximum_exo_lag+maximum_exo_lead+1,1]), params0, ys, 1);
        %g1 is [endo_nbr by (dynamicvar_nbr+exo_nbr)] first derivative (wrt all endogenous, exogenous and auxiliary variables) of dynamic model equations, i.e. df/d[yy0;ex0], in DR order
        %g2 is [endo_nbr by (dynamicvar_nbr+exo_nbr)^2] second derivative (wrt all endogenous, exogenous and auxiliary variables) of dynamic model equations, i.e. d(df/d[yy0;ex0])/d[yy0;ex0], in DR order
    end
    yy0ex0_nbr = sqrt(size(g2,2)); % number of dynamic variables + exogenous variables (length(yy0)+length(ex0))
    dYss = dys(oo.dr.order_var, indpmodel); %focus only on model parameters, note dys is in declaration order, dYss is in DR-order
    dyy0 = dys(I,:);
    yy0_nbr = max(max(lead_lag_incidence)); % retrieve the number of states excluding columns for shocks
    % Computation of dg1, i.e. first derivative (wrt. parameters) of Jacobian (wrt endogenous and auxilary variables) of dynamic model using the implicit function theorem
    %   Let g1 denote the Jacobian of dynamic model equations, i.e. g1 = df/d[yy0ex0], evaluated at the steady state
    %   Let dg1 denote the first-derivative (wrt parameters) of g1 evaluated at the steady state
    %   Note that g1 is a function of both the parameters and of the steady state, which also depends on the parameters.
    %   Hence, implicitly g1=g1(p,yy0ex0(p)) and dg1 consists of two parts (see Ratto and Iskrev (2011) formula 7):
    %   (1) direct derivative wrt to parameters given by the preprocessor, i.e. gp
    %   and
    %   (2) contribution of derivative of steady state (wrt parameters), i.e. g2*dyy0
    %   Note that in a stochastic context ex0 is always zero and hence can be skipped in the computations
    dg1_part2 = gp*0; %initialize part 2, it has dimension [endo_nbr by (dynamicvar_nbr+exo_nbr) by param_nbr]
    for j = 1:endo_nbr
        [II, JJ] = ind2sub([yy0ex0_nbr yy0ex0_nbr], find(g2(j,:)));
        %g2 is [endo_nbr by (dynamicvar_nbr+exo_nbr)^2]
        for i = 1:yy0ex0_nbr
            is = find(II==i);
            is = is(find(JJ(is)<=yy0_nbr));  %focus only on yy0 derivatives as ex0 variables are 0 in a stochastic context
            if ~isempty(is)
                tmp_g2 = full(g2(j,find(g2(j,:))));
                dg1_part2(j,i,:) = tmp_g2(is)*dyy0(JJ(is),:); %put into tensor notation
            end
        end
    end
    dg1 = gp + dg1_part2;     %dg is sum of two parts due to implicit function theorem
    dg1 = dg1(:,:,indpmodel); %focus only on model parameters

    if nargout > 6
        % Computation of d2g1, i.e. second derivative (wrt. parameters) of Jacobian (wrt endogenous and auxilary variables) of dynamic model using the implicit function theorem
        %   Let g1 denote the Jacobian of dynamic model equations, i.e. g1 = df/d[yy0ex0], evaluated at the steady state
        %   Let d2g1 denote the second-derivative (wrt parameters) of g1
        %   Note that g1 is a function of both the parameters and of the steady state, which also depends on the parameters.
        %   Hence, implicitly g1=g1(p,yy0ex0(p)) and the first derivative is given by dg1 = gp + g2*dyy0ex0 (see above)
        % Accordingly, d2g1, the second-derivative (wrt parameters), consists of five parts (ignoring transposes, see Ratto and Iskrev (2011) formula 16)
        %   (1) d(gp)/dp                              = gpp
        %   (2) d(gp)/dyy0ex0*d(yy0ex0)/dp            = hp  * dyy0ex0
        %   (3) d(g2)/dp * dyy0ex0                    = hp  * dyy0ex0
        %   (4) d(g2)/dyy0ex0*d(dyy0ex0)/dp * dyy0ex0 = g3  * dyy0ex0  * dyy0ex0
        %   (5) g2 * d(dyy0ex0)/dp                    = g2  * d2yy0ex0
        %   Note that part 2 and 3 are equivalent besides the use of transpose (see Ratto and Iskrev (2011) formula 16)
        d2g1_full = sparse(endo_nbr*yy0ex0_nbr, param_nbr*param_nbr); %initialize
        dyy0ex0 = sparse([dyy0; zeros(yy0ex0_nbr-yy0_nbr,param_nbr)]);     %Jacobian (wrt model parameters) of steady state of dynamic (endogenous and auxiliary) and exogenous variables
        
        g3_tmp = reshape(g3,[endo_nbr*yy0ex0_nbr*yy0ex0_nbr yy0ex0_nbr]); 
        d2g1_part4_left = sparse(endo_nbr*yy0ex0_nbr*yy0ex0_nbr,param_nbr);
        for j = 1:param_nbr
            %compute first two terms of part 4
            d2g1_part4_left(:,j) = g3_tmp*dyy0ex0(:,j);
        end
       
        for j=1:endo_nbr
            %Note that in the following we focus only on dynamic variables as exogenous variables are 0 by construction in a stochastic setting
            d2g1_part5 = reshape(g2(j,:), [yy0ex0_nbr yy0ex0_nbr]);
            d2g1_part5 = d2g1_part5(:,1:yy0_nbr)*reshape(d2ys(I,:,:),[yy0_nbr,param_nbr*param_nbr]);
            for i=1:yy0ex0_nbr
                ind_part4 = sub2ind([endo_nbr yy0ex0_nbr yy0ex0_nbr], ones(yy0ex0_nbr,1)*j ,ones(yy0ex0_nbr,1)*i, (1:yy0ex0_nbr)');
                d2g1_part4 = (d2g1_part4_left(ind_part4,:))'*dyy0ex0;
                d2g1_part2_and_part3 = (get_hess_deriv(hp,j,i,yy0ex0_nbr,param_nbr))'*dyy0ex0;
                d2g1_part1 = get_2nd_deriv_mat(gpp,j,i,param_nbr);
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
        clear d2g1_full;
    end
end
% clear variables that are not used any more
clear rp_static g1_static 
clear ys dys dyy0 dyy0ex0 
clear dg1_part2 tmp_g2
clear g2 gp rpp_static g2_static gp_static d2ys
clear hp g3 g3_tmp gpp
clear ind_d2g1 ind_d2g1_tmp ind_part4 i j i1 i2 ig1 ig2 I II JJ ip1 ip2 is

% Construct nonzero derivatives wrt to t+1, t, and t-1 variables using kstate
klen = maximum_endo_lag + maximum_endo_lead + 1; %total length
k11 = lead_lag_incidence(find([1:klen] ~= maximum_endo_lag+1),:);
g1nonzero = g1(:,nonzeros(k11'));
dg1nonzero = dg1(:,nonzeros(k11'),:);
if nargout > 6
    indind = ismember(d2g1(:,2),nonzeros(k11'));
    tmp = d2g1(indind,:);
    d2g1nonzero = tmp;
    for j = 1:size(tmp,1)
        inxinx = find(nonzeros(k11')==tmp(j,2));
        d2g1nonzero(j,2) = inxinx;
    end
end
kstate = oo.dr.kstate;

% Construct nonzero derivatives wrt to t+1, i.e. GAM1=-f_{y^+} in Villemot (2011)
GAM1 = zeros(endo_nbr,endo_nbr);
dGAM1 = zeros(endo_nbr,endo_nbr,modparam_nbr);
k1 = find(kstate(:,2) == maximum_endo_lag+2 & kstate(:,3));
GAM1(:, kstate(k1,1)) = -g1nonzero(:,kstate(k1,3));
dGAM1(:, kstate(k1,1), :) = -dg1nonzero(:,kstate(k1,3),:);
if nargout > 6
    indind = ismember(d2g1nonzero(:,2),kstate(k1,3));
    tmp = d2g1nonzero(indind,:);
    tmp(:,end)=-tmp(:,end);
    d2GAM1 = tmp;
    for j = 1:size(tmp,1)
        inxinx = (kstate(k1,3)==tmp(j,2));
        d2GAM1(j,2) = kstate(k1(inxinx),1);
    end
end

% Construct nonzero derivatives wrt to t, i.e. GAM0=f_{y^0} in Villemot (2011)
[~,cols_b,cols_j] = find(lead_lag_incidence(maximum_endo_lag+1, oo.dr.order_var));
GAM0 = zeros(endo_nbr,endo_nbr);
dGAM0 = zeros(endo_nbr,endo_nbr,modparam_nbr);
GAM0(:,cols_b) = g1(:,cols_j);
dGAM0(:,cols_b,:) = dg1(:,cols_j,:);
if nargout > 6
    indind = ismember(d2g1(:,2),cols_j);
    tmp = d2g1(indind,:);
    d2GAM0 = tmp;
    for j = 1:size(tmp,1)
        inxinx = (cols_j==tmp(j,2));
        d2GAM0(j,2) = cols_b(inxinx);
    end
end

% Construct nonzero derivatives wrt to t-1, i.e. GAM2=-f_{y^-} in Villemot (2011)
k2 = find(kstate(:,2) == maximum_endo_lag+1 & kstate(:,4));
GAM2 = zeros(endo_nbr,endo_nbr);
dGAM2 = zeros(endo_nbr,endo_nbr,modparam_nbr);
GAM2(:, kstate(k2,1)) = -g1nonzero(:,kstate(k2,4));
dGAM2(:, kstate(k2,1), :) = -dg1nonzero(:,kstate(k2,4),:);
if nargout > 6
    indind = ismember(d2g1nonzero(:,2),kstate(k2,4));
    tmp = d2g1nonzero(indind,:);
    tmp(:,end) = -tmp(:,end);
    d2GAM2 = tmp;
    for j = 1:size(tmp,1)
        inxinx = (kstate(k2,4)==tmp(j,2));
        d2GAM2(j,2) = kstate(k2(inxinx),1);
    end
end

% Construct nonzero derivatives wrt to u_t, i.e. GAM3=-f_{u} in Villemot (2011)
GAM3 = -g1(:,length(yy0)+1:end);
dGAM3 = -dg1(:,length(yy0)+1:end,:);
if nargout > 6
    cols_ex = [length(yy0)+1:size(g1,2)];
    indind = ismember(d2g1(:,2),cols_ex);
    tmp = d2g1(indind,:);
    tmp(:,end) = -tmp(:,end);
    d2GAM3 = tmp;
    for j = 1:size(tmp,1)
        inxinx = find(cols_ex==tmp(j,2));
        d2GAM3(j,2) = inxinx;
    end
    clear d2g1 d2g1nonzero tmp
end
clear cols_b cols_ex cols_j k1 k11 k2 klen kstate
clear g1nonzero dg1nonzero g1 yy0

%% Construct first derivative of Sigma_e
dSigma_e = zeros(exo_nbr,exo_nbr,totparam_nbr); %initialize
% note that derivatives wrt model parameters are zero by construction
% Compute first derivative of Sigma_e wrt stderr parameters (these come first)
if ~isempty(indpstderr)
    for jp = 1:stderrparam_nbr
        dSigma_e(indpstderr(jp),indpstderr(jp),jp) = 2*stderr_e0(indpstderr(jp));
        if isdiag(Sigma_e0) == 0 % if there are correlated errors add cross derivatives
            indotherex0 = 1:exo_nbr;
            indotherex0(indpstderr(jp)) = [];
            for kk = indotherex0
                dSigma_e(indpstderr(jp), kk, jp) = Corr_e0(indpstderr(jp),kk)*stderr_e0(kk);
                dSigma_e(kk, indpstderr(jp), jp) = dSigma_e(indpstderr(jp), kk, jp); %symmetry
            end
        end
    end
end
% Compute first derivative of Sigma_e wrt corr parameters (these come second)
if ~isempty(indpcorr)
    for jp = 1:corrparam_nbr
        dSigma_e(indpcorr(jp,1),indpcorr(jp,2),stderrparam_nbr+jp) = stderr_e0(indpcorr(jp,1))*stderr_e0(indpcorr(jp,2));
        dSigma_e(indpcorr(jp,2),indpcorr(jp,1),stderrparam_nbr+jp) = dSigma_e(indpcorr(jp,1),indpcorr(jp,2),stderrparam_nbr+jp); %symmetry
    end
end

%% Construct second derivative of Sigma_e
if nargout > 6
    % note that derivatives wrt (mod x mod) and (corr x corr) parameters
    % are zero by construction; hence we only need to focus on (stderr x stderr), and (stderr x corr)
    d2Sigma_e = zeros(exo_nbr,exo_nbr,totparam_nbr^2); %initialize full matrix, even though we'll reduce it later on to unique upper triangular values
    % Compute upper triangular values of Hessian of Sigma_e wrt (stderr x stderr) parameters
    if ~isempty(indp2stderrstderr)
        for jp = 1:stderrparam_nbr
            for ip = 1:jp                
                if jp == ip %same stderr parameters
                    d2Sigma_e(indpstderr(jp),indpstderr(jp),indp2stderrstderr(ip,jp)) = 2;
                else %different stderr parameters
                    if isdiag(Sigma_e0) == 0 % if there are correlated errors
                        d2Sigma_e(indpstderr(jp),indpstderr(ip),indp2stderrstderr(ip,jp)) = Corr_e0(indpstderr(jp),indpstderr(ip));
                        d2Sigma_e(indpstderr(ip),indpstderr(jp),indp2stderrstderr(ip,jp)) = Corr_e0(indpstderr(jp),indpstderr(ip)); %symmetry
                    end
                end
            end
        end
    end
    % Compute upper triangular values of Hessian of Sigma_e wrt (stderr x corr) parameters
    if ~isempty(indp2stderrcorr)
        for jp = 1:stderrparam_nbr
            for ip = 1:corrparam_nbr                
                if indpstderr(jp) == indpcorr(ip,1) %if stderr equal to first index of corr parameter, derivative is equal to stderr corresponding to second index
                    d2Sigma_e(indpstderr(jp),indpcorr(ip,2),indp2stderrcorr(jp,ip)) = stderr_e0(indpcorr(ip,2));
                    d2Sigma_e(indpcorr(ip,2),indpstderr(jp),indp2stderrcorr(jp,ip)) = stderr_e0(indpcorr(ip,2)); % symmetry
                end
                if indpstderr(jp) == indpcorr(ip,2) %if stderr equal to second index of corr parameter, derivative is equal to stderr corresponding to first index
                    d2Sigma_e(indpstderr(jp),indpcorr(ip,1),indp2stderrcorr(jp,ip)) = stderr_e0(indpcorr(ip,1));
                    d2Sigma_e(indpcorr(ip,1),indpstderr(jp),indp2stderrcorr(jp,ip)) = stderr_e0(indpcorr(ip,1)); % symmetry
                end
            end
        end
    end
    d2Sigma_e = d2Sigma_e(:,:,indp2tottot2); %focus on upper triangular hessian values
end


if kronflag == 1
    % The following derivations are based on Iskrev (2010) and its online appendix A. 
    % Basic idea is to make use of the implicit function theorem.
    % Let F = GAM0*A - GAM1*A*A - GAM2 = 0
    % Note that F is a function of parameters p and A, which is also a
    % function of p,therefore, F = F(p,A(p)), and hence, 
    % dF = Fp + dF_dA*dA or dA = - Fp/dF_dA
    
    % Some auxiliary matrices
    I_endo = speye(endo_nbr);
    I_exo = speye(exo_nbr);

    % Reshape to write derivatives in the Magnus and Neudecker style, i.e. dvec(X)/dp
    dGAM0 = reshape(dGAM0, endo_nbr^2, modparam_nbr);
    dGAM1 = reshape(dGAM1, endo_nbr^2, modparam_nbr);
    dGAM2 = reshape(dGAM2, endo_nbr^2, modparam_nbr);
    dGAM3 = reshape(dGAM3, endo_nbr*exo_nbr, modparam_nbr);
    dSigma_e = reshape(dSigma_e, exo_nbr^2, totparam_nbr);
    
    % Compute dA via implicit function
    dF_dA = kron(I_endo,GAM0) - kron(A',GAM1) - kron(I_endo,GAM1*A); %equation 31 in Appendix A of Iskrev (2010)
    Fp = kron(A',I_endo)*dGAM0 - kron( (A')^2,I_endo)*dGAM1 - dGAM2; %equation 32 in Appendix A of Iskrev (2010)
    dA = -dF_dA\Fp;
    
    % Compute dB from expressions 33 in Iskrev (2010) Appendix A
    MM = GAM0-GAM1*A; %this corresponds to matrix M in Ratto and Iskrev (2011, page 6) and will be used if nargout > 6 below
    invMM = MM\eye(endo_nbr);
    dB = - kron( (invMM*GAM3)' , invMM ) * ( dGAM0 - kron( A' , I_endo ) * dGAM1 - kron( I_endo , GAM1 ) * dA ) + kron( I_exo, invMM ) * dGAM3 ;
    dBt = commutation(endo_nbr, exo_nbr)*dB; %transose of derivative using the commutation matrix
    
    % Add derivatives for stderr and corr parameters, which are zero by construction
    dA = [zeros(endo_nbr^2, stderrparam_nbr+corrparam_nbr) dA];
    dB = [zeros(endo_nbr*exo_nbr, stderrparam_nbr+corrparam_nbr) dB];    
    dBt = [zeros(endo_nbr*exo_nbr, stderrparam_nbr+corrparam_nbr) dBt];
    
    % Compute dOm = dvec(B*Sig*B') from expressions 34 in Iskrev (2010) Appendix A    
    dOm = kron(I_endo,B*Sigma_e0)*dBt + kron(B,B)*dSigma_e + kron(B*Sigma_e0,I_endo)*dB;
    
    % Put into tensor notation
    dA   = reshape(dA,   endo_nbr, endo_nbr, totparam_nbr);
    dB   = reshape(dB,   endo_nbr, exo_nbr,  totparam_nbr);
    dOm  = reshape(dOm,  endo_nbr, endo_nbr, totparam_nbr);
    dSigma_e = reshape(dSigma_e, exo_nbr,  exo_nbr,  totparam_nbr);
    if nargout > 6
        % Put back into tensor notation as these will be reused later
        dGAM0 = reshape(dGAM0, endo_nbr, endo_nbr, modparam_nbr);
        dGAM1 = reshape(dGAM1, endo_nbr, endo_nbr, modparam_nbr);
        dGAM2 = reshape(dGAM2, endo_nbr, endo_nbr, modparam_nbr);
        dGAM3 = reshape(dGAM3, endo_nbr, exo_nbr,  modparam_nbr);        
        dAA   = dA(:, :, stderrparam_nbr+corrparam_nbr+1:end); %this corresponds to matrix dA in Ratto and Iskrev (2011, page 6), i.e. derivative of A with respect to model parameters only in tensor notation
        dBB   = dB(:, :, stderrparam_nbr+corrparam_nbr+1:end); %dBB is for all endogenous variables, whereas dB is only for selected variables
        N = -GAM1; %this corresponds to matrix N in Ratto and Iskrev (2011, page 6)
        P = A;     %this corresponds to matrix P in Ratto and Iskrev (2011, page 6)
    end
    
    % Focus only on selected variables
    dYss = dYss(indvar,:);
    dA   = dA(indvar,indvar,:);
    dB   = dB(indvar,:,:);
    dOm  = dOm(indvar,indvar,:);
    
elseif (kronflag == 0 || kronflag == -2)
    % generalized sylvester equation solves MM*dAA+N*dAA*P=Q from Ratto and Iskrev (2011) equation 11 where
    % dAA is derivative of A with respect to model parameters only in tensor notation
    MM = (GAM0-GAM1*A);
    N = -GAM1;
    P = A;
    Q_rightpart = zeros(endo_nbr,endo_nbr,modparam_nbr); %initialize
    Q = Q_rightpart; %initialize and compute matrix Q in Ratto and Iskrev (2011, page 6)
    for j = 1:modparam_nbr
        Q_rightpart(:,:,j) = (dGAM0(:,:,j)-dGAM1(:,:,j)*A);
        Q(:,:,j) = dGAM2(:,:,j)-Q_rightpart(:,:,j)*A;
    end
    %use iterated generalized sylvester equation to compute dAA
    dAA = sylvester3(MM,N,P,Q);
    flag = 1; icount = 0;
    while flag && icount < 4
        [dAA, flag] = sylvester3a(dAA,MM,N,P,Q);
        icount = icount+1;
    end
    
    %stderr parameters come first, then corr parameters, model parameters come last
    %note that stderr and corr derivatives are:
    % - zero by construction for A and B
    % - depend only on dSig for Om
    dOm = zeros(var_nbr,  var_nbr, totparam_nbr);
    dA  = zeros(var_nbr,  var_nbr, totparam_nbr);
    dB  = zeros(var_nbr, exo_nbr, totparam_nbr);
    if nargout > 6
        dBB  = zeros(endo_nbr, exo_nbr, modparam_nbr); %dBB is always for all endogenous variables, whereas dB is only for selected variables
    end
    
    %compute derivative of Om=B*Sig*B' that depends on Sig (other part is added later)
    if ~isempty(indpstderr)
        for j = 1:stderrparam_nbr
            BSigjBt = B*dSigma_e(:,:,j)*B';            
            dOm(:,:,j) = BSigjBt(indvar,indvar);            
        end
    end
    if ~isempty(indpcorr)
        for j = 1:corrparam_nbr
            BSigjBt = B*dSigma_e(:,:,stderrparam_nbr+j)*B';
            dOm(:,:,stderrparam_nbr+j) = BSigjBt(indvar,indvar);
        end
    end
    
    %compute derivative of B and the part of Om=B*Sig*B' that depends on B (other part is computed above)
    invMM = inv(MM);
    for j = 1:modparam_nbr        
        dAAj = dAA(:,:,j);
        dBj = invMM * ( dGAM3(:,:,j) - (Q_rightpart(:,:,j) -GAM1*dAAj ) * B ); %equation 14 in Ratto and Iskrev (2011), except in the paper there is a typo as the last B is missing
        dOmj = dBj*Sigma_e0*B'+B*Sigma_e0*dBj';
        %store derivatives in tensor notation
        dA(:, :, stderrparam_nbr+corrparam_nbr+j) = dAAj(indvar,indvar);
        dB(:, :, stderrparam_nbr+corrparam_nbr+j) = dBj(indvar,:);
        dOm(:, :, stderrparam_nbr+corrparam_nbr+j) = dOmj(indvar,indvar);
        if nargout > 6
            dBB(:, :, j) = dBj;
        end
    end    
    dYss = dYss(indvar,:); % Focus only on relevant variables
end

%% Compute second-order derivatives (wrt params) of solution matrices using generalized sylvester equations, see equations 17 and 18 in Ratto and Iskrev (2011)
if nargout > 6
    % solves MM*d2AA+N*d2AA*P = QQ where d2AA are second order derivatives (wrt model parameters) of A
    d2Yss = d2Yss(indvar,:,:);
    QQ = zeros(endo_nbr,endo_nbr,floor(sqrt(modparam_nbr2)));
    jcount=0;
    cumjcount=0;
    jinx = [];
    x2x=sparse(endo_nbr*endo_nbr,modparam_nbr2);
    for i=1:modparam_nbr
        for j=1:i
            elem1 = (get_2nd_deriv(d2GAM0,endo_nbr,endo_nbr,j,i)-get_2nd_deriv(d2GAM1,endo_nbr,endo_nbr,j,i)*A);
            elem1 = get_2nd_deriv(d2GAM2,endo_nbr,endo_nbr,j,i)-elem1*A;
            elemj0 = dGAM0(:,:,j)-dGAM1(:,:,j)*A;
            elemi0 = dGAM0(:,:,i)-dGAM1(:,:,i)*A;
            elem2 = -elemj0*dAA(:,:,i)-elemi0*dAA(:,:,j);
            elem2 = elem2 + ( dGAM1(:,:,j)*dAA(:,:,i) + dGAM1(:,:,i)*dAA(:,:,j) )*A;
            elem2 = elem2 + GAM1*( dAA(:,:,i)*dAA(:,:,j) + dAA(:,:,j)*dAA(:,:,i));
            jcount=jcount+1;
            jinx = [jinx; [j i]];
            QQ(:,:,jcount) = elem1+elem2;
            if jcount==floor(sqrt(modparam_nbr2)) || (j*i)==modparam_nbr^2
                if (j*i)==modparam_nbr^2
                    QQ = QQ(:,:,1:jcount);
                end
                xx2=sylvester3(MM,N,P,QQ);
                flag=1;
                icount=0;
                while flag && icount<4
                    [xx2, flag]=sylvester3a(xx2,MM,N,P,QQ);
                    icount = icount + 1;
                end
                x2x(:,cumjcount+1:cumjcount+jcount)=reshape(xx2,[endo_nbr*endo_nbr jcount]);
                cumjcount=cumjcount+jcount;
                jcount = 0;
                jinx = [];
            end
        end
    end
    clear d xx2;
    jcount = 0;
    icount = 0;
    cumjcount = 0;
    MAX_DIM_MAT = 100000000;
    ncol = max(1,floor(MAX_DIM_MAT/(8*var_nbr*(var_nbr+1)/2)));
    ncol = min(ncol, totparam_nbr2);
    d2A = sparse(var_nbr*var_nbr,totparam_nbr2);
    d2Om = sparse(var_nbr*(var_nbr+1)/2,totparam_nbr2);
    d2A_tmp = zeros(var_nbr*var_nbr,ncol);
    d2Om_tmp = zeros(var_nbr*(var_nbr+1)/2,ncol);
    tmpDir = CheckPath('tmp_derivs',dname);
    offset = stderrparam_nbr+corrparam_nbr;
    %     d2B = zeros(m,n,tot_param_nbr,tot_param_nbr);    
    for j=1:totparam_nbr
        for i=1:j
            jcount=jcount+1;
            if j<=offset %stderr and corr parameters                
                    y = B*d2Sigma_e(:,:,jcount)*B';
                    d2Om_tmp(:,jcount) = dyn_vech(y(indvar,indvar));
            else %model parameters
                jind = j-offset;
                iind = i-offset;
                if i<=offset
                    y = dBB(:,:,jind)*dSigma_e(:,:,i)*B'+B*dSigma_e(:,:,i)*dBB(:,:,jind)';
                    %                     y(abs(y)<1.e-8)=0;
                    d2Om_tmp(:,jcount) = dyn_vech(y(indvar,indvar));
                else
                    icount=icount+1;
                    dAAj = reshape(x2x(:,icount),[endo_nbr endo_nbr]);
                    %                     x = get_2nd_deriv(x2x,m,m,iind,jind);%xx2(:,:,jcount);
                    elem1 = (get_2nd_deriv(d2GAM0,endo_nbr,endo_nbr,iind,jind)-get_2nd_deriv(d2GAM1,endo_nbr,endo_nbr,iind,jind)*A);
                    elem1 = elem1 -( dGAM1(:,:,jind)*dAA(:,:,iind) + dGAM1(:,:,iind)*dAA(:,:,jind) );
                    elemj0 = dGAM0(:,:,jind)-dGAM1(:,:,jind)*A-GAM1*dAA(:,:,jind);
                    elemi0 = dGAM0(:,:,iind)-dGAM1(:,:,iind)*A-GAM1*dAA(:,:,iind);
                    elem0 = elemj0*dBB(:,:,iind)+elemi0*dBB(:,:,jind);
                    y = invMM * (get_2nd_deriv(d2GAM3,endo_nbr,exo_nbr,iind,jind)-elem0-(elem1-GAM1*dAAj)*B);
                    %         d2B(:,:,j+length(indexo),i+length(indexo)) = y;
                    %         d2B(:,:,i+length(indexo),j+length(indexo)) = y;
                    y = y*Sigma_e0*B'+B*Sigma_e0*y'+ ...
                        dBB(:,:,jind)*Sigma_e0*dBB(:,:,iind)'+dBB(:,:,iind)*Sigma_e0*dBB(:,:,jind)';
                    %                     x(abs(x)<1.e-8)=0;
                    d2A_tmp(:,jcount) = vec(dAAj(indvar,indvar));
                    %                     y(abs(y)<1.e-8)=0;
                    d2Om_tmp(:,jcount) = dyn_vech(y(indvar,indvar));
                end
            end
            if jcount==ncol || i*j==totparam_nbr^2
                d2A(:,cumjcount+1:cumjcount+jcount) = d2A_tmp(:,1:jcount);
                %         d2A(:,:,j+length(indexo),i+length(indexo)) = x;
                %         d2A(:,:,i+length(indexo),j+length(indexo)) = x;
                d2Om(:,cumjcount+1:cumjcount+jcount) = d2Om_tmp(:,1:jcount);
                %         d2Om(:,:,j+length(indexo),i+length(indexo)) = y;
                %         d2Om(:,:,i+length(indexo),j+length(indexo)) = y;
                save([tmpDir filesep 'd2A_' int2str(cumjcount+1) '_' int2str(cumjcount+jcount) '.mat'],'d2A')
                save([tmpDir filesep 'd2Om_' int2str(cumjcount+1) '_'  int2str(cumjcount+jcount) '.mat'],'d2Om')
                cumjcount = cumjcount+jcount;
                jcount=0;
                %         d2A = sparse(m1*m1,tot_param_nbr*(tot_param_nbr+1)/2);
                %         d2Om = sparse(m1*(m1+1)/2,tot_param_nbr*(tot_param_nbr+1)/2);
                d2A_tmp = zeros(var_nbr*var_nbr,ncol);
                d2Om_tmp = zeros(var_nbr*(var_nbr+1)/2,ncol);
            end
        end
    end
end

return

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
% - npar: scalar                                     Number of model parameters, i.e. equals M_.param_nbr
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

function g22 = get_all_2nd_derivs(gpp,m,n,npar,fsparse)

if nargin==4 || isempty(fsparse)
    fsparse=0;
end
if fsparse
    g22=sparse(m*n,npar*npar);
else
    g22=zeros(m,n,npar,npar);
end
% c=ones(npar,npar);
% c=triu(c);
% ic=find(c);

for is=1:length(gpp)
    %     d=zeros(npar,npar);
    %     d(gpp(is,3),gpp(is,4))=1;
    %     indx = find(ic==find(d));
    if fsparse
        g22(sub2ind([m,n],gpp(is,1),gpp(is,2)),sub2ind([npar,npar],gpp(is,3),gpp(is,4)))=gpp(is,5);
    else
        g22(gpp(is,1),gpp(is,2),gpp(is,3),gpp(is,4))=gpp(is,5);
    end
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

for is=1:length(rpp)
    % Keep symmetry in hessian, hence 2 and 3 as well as 3 and 2, i.e. d2f/(dp1 dp2) = d2f/(dp2 dp1)
    r22(rpp(is,1),rpp(is,2),rpp(is,3))=rpp(is,4);
    r22(rpp(is,1),rpp(is,3),rpp(is,2))=rpp(is,4);
end

return

function h2 = get_all_hess_derivs(hp,r,m,npar)

h2=zeros(r,m,m,npar);

for is=1:length(hp)
    h2(hp(is,1),hp(is,2),hp(is,3),hp(is,4))=hp(is,5);
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
% - m:    scalar                                     Number of dynamic model variables + exogenous vars, i.e. dynamicvar_nbr + exo_nbr
% - npar: scalar                                     Number of model parameters, i.e. equals M_.param_nbr
%
% output:
% h2: [(dynamicvar_nbr + exo_nbr) by M_.param_nbr] Jacobian matrix (wrt parameters) of dynamic Hessian
%                                                              rows: second dynamic or exogenous variables in Hessian of specific model equation of the dynamic model
%                                                              columns: parameters

h2=zeros(m,npar);
is1=find(hp(:,1)==i);
is=is1(find(hp(is1,2)==j));

if ~isempty(is)
    h2(sub2ind([m,npar],hp(is,3),hp(is,4)))=hp(is,5)';
end

return
