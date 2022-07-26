function [alphahat,etahat,epsilonhat,ahat,SteadyState,trend_coeff,aK,T,R,P,PK,decomp,trend_addition,state_uncertainty,M_,oo_,bayestopt_] = DsgeSmoother(xparam1,gend,Y,data_index,missing_value,M_,oo_,options_,bayestopt_,estim_params_,varargin)
% Estimation of the smoothed variables and innovations.
%
% INPUTS
%   o xparam1       [double]   (p*1) vector of (estimated) parameters.
%   o gend          [integer]  scalar specifying the number of observations ==> varargin{1}.
%   o Y             [double]   (n*T) matrix of data.
%   o data_index    [cell]      1*smpl cell of column vectors of indices.
%   o missing_value 1 if missing values, 0 otherwise
%   o M_            [structure] decribing the model
%   o oo_           [structure] storing the results
%   o options_      [structure] describing the options
%   o bayestopt_    [structure] describing the priors
%   o estim_params_ [structure] characterizing parameters to be estimated
%
% OUTPUTS
%   o alphahat      [double]  (m*T) matrix, smoothed endogenous variables (a_{t|T})  (decision-rule order)
%   o etahat        [double]  (r*T) matrix, smoothed structural shocks (r>=n is the number of shocks).
%   o epsilonhat    [double]  (n*T) matrix, smoothed measurement errors.
%   o ahat          [double]  (m*T) matrix, updated (endogenous) variables (a_{t|t}) (decision-rule order)
%   o SteadyState   [double]  (m*1) vector specifying the steady state level of each endogenous variable (declaration order)
%   o trend_coeff   [double]  (n*1) vector, parameters specifying the slope of the trend associated to each observed variable.
%   o aK            [double]  (K,n,T+K) array, k (k=1,...,K) steps ahead
%                                   filtered (endogenous) variables  (decision-rule order)
%   o T and R       [double]  Matrices defining the state equation (T is the (m*m) transition matrix).
%   o P:            (m*m*(T+1)) 3D array of one-step ahead forecast error variance
%                       matrices (decision-rule order)
%   o PK:           (K*m*m*(T+K)) 4D array of k-step ahead forecast error variance
%                       matrices (meaningless for periods 1:d) (decision-rule order)
%   o decomp        (K*m*r*(T+K)) 4D array of shock decomposition of k-step ahead
%                       filtered variables (decision-rule order)
%   o trend_addition [double] (n*T) pure trend component; stored in options_.varobs order
%   o state_uncertainty [double] (K,K,T) array, storing the uncertainty
%                                   about the smoothed state (decision-rule order)
%   o M_            [structure] decribing the model
%   o oo_           [structure] storing the results
%   o bayestopt_    [structure] describing the priors
%
% Notes:
%   m:  number of endogenous variables (M_.endo_nbr)
%   T:  number of Time periods (options_.nobs)
%   r:  number of strucural shocks (M_.exo_nbr)
%   n:  number of observables (length(options_.varobs))
%   K:  maximum forecast horizon (max(options_.nk))
%
%   To get variables that are stored in decision rule order in order of declaration
%   as in M_.endo_names, ones needs code along the lines of:
%   variables_declaration_order(dr.order_var,:) = alphahat
%
%   Defines bayestopt_.mf = bayestopt_.smoother_mf (positions of observed variables
%   and requested smoothed variables in decision rules (decision rule order)) and
%   passes it back via global variable
%
% ALGORITHM
%   Diffuse Kalman filter (Durbin and Koopman)
%
% SPECIAL REQUIREMENTS
%   None

% Copyright (C) 2006-2020 Dynare Team
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

alphahat        = [];
etahat  = [];
epsilonhat      = [];
ahat          = [];
SteadyState   = [];
trend_coeff   = [];
aK            = [];
T             = [];
R             = [];
P             = [];
PK            = [];
decomp        = [];
vobs            = length(options_.varobs);
smpl          = size(Y,2);

if ~isempty(xparam1) %not calibrated model
    M_ = set_all_parameters(xparam1,estim_params_,M_);
end

%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------
length_varargin=length(varargin);
if ~options_.smoother_redux
    
    %store old setting of restricted var_list
    oldoo.restrict_var_list = oo_.dr.restrict_var_list;
    oldoo.restrict_columns = oo_.dr.restrict_columns;
    oo_.dr.restrict_var_list = bayestopt_.smoother_var_list;
    oo_.dr.restrict_columns = bayestopt_.smoother_restrict_columns;
    
    [T,R,SteadyState,info,M_,oo_] = dynare_resolve(M_,options_,oo_);
    
    %get location of observed variables and requested smoothed variables in
    %decision rules
    bayestopt_.mf = bayestopt_.smoother_var_list(bayestopt_.smoother_mf);
    
else
    if ~options_.occbin.smoother.status
        [T,R,SteadyState,info,M_,oo_] = dynare_resolve(M_,options_,oo_,'restrict');
    else
        [T,R,SteadyState,info,M_,oo_,~,~,~, T0, R0] = ...
            occbin.dynare_resolve(M_,options_,oo_,[],'restrict');
        varargin{length_varargin+1}=T0;
        varargin{length_varargin+2}=R0;
    end
    bayestopt_.mf = bayestopt_.mf1;
end
if options_.occbin.smoother.status
    occbin_info.status = true;
    occbin_info.info= [{options_,oo_,M_} varargin];
else
    occbin_info.status = false;    
end

if info~=0
    print_info(info,options_.noprint, options_);
    return
end

if options_.noconstant
    constant = zeros(vobs,1);
else
    if options_.loglinear
        constant = log(SteadyState(bayestopt_.mfys));
    else
        constant = SteadyState(bayestopt_.mfys);
    end
end
trend_coeff = zeros(vobs,1);
if bayestopt_.with_trend == 1
    [trend_addition, trend_coeff] =compute_trend_coefficients(M_,options_,vobs,gend);
    trend = constant*ones(1,gend)+trend_addition;
else
    trend_addition=zeros(size(constant,1),gend);
    trend = constant*ones(1,gend);
end
start = options_.presample+1;
np    = size(T,1);
mf    = bayestopt_.mf;
% ------------------------------------------------------------------------------
%  3. Initial condition of the Kalman filter
% ------------------------------------------------------------------------------
%
%  Here, Pinf and Pstar are determined. If the model is stationary, determine
%  Pstar as the solution of the Lyapunov equation and set Pinf=[] (Notation follows
%  Koopman/Durbin (2003), Journal of Time Series Analysis 24(1))
%
Q = M_.Sigma_e;
H = M_.H;

if isequal(H,0)
    H = zeros(vobs,vobs);
end

Z = zeros(vobs,size(T,2));
for i=1:vobs
    Z(i,mf(i)) = 1;
end

expanded_state_vector_for_univariate_filter=0;
kalman_algo = options_.kalman_algo;
if options_.lik_init == 1               % Kalman filter
    if kalman_algo ~= 2
        kalman_algo = 1;
    end
    Pstar=lyapunov_solver(T,R,Q,options_);
    Pinf        = [];
elseif options_.lik_init == 2           % Old Diffuse Kalman filter
    if kalman_algo ~= 2
        kalman_algo = 1;
    end
    Pstar = options_.Harvey_scale_factor*eye(np);
    Pinf        = [];
elseif options_.lik_init == 3           % Diffuse Kalman filter
    if kalman_algo ~= 4
        kalman_algo = 3;
    else
        if ~all(all(abs(H-diag(diag(H)))<1e-14))% ie, the covariance matrix is not diagonal...
            %Augment state vector (follows Section 6.4.3 of DK (2012))
            expanded_state_vector_for_univariate_filter=1;
            T  = blkdiag(T,zeros(vobs));
            np    = size(T,1);
            Q   = blkdiag(Q,H);
            R  = blkdiag(R,eye(vobs));
            H   = zeros(vobs,vobs);
            Z   = [Z, eye(vobs)];
        end
    end
    [Pstar,Pinf] = compute_Pinf_Pstar(mf,T,R,Q,options_.qz_criterium);
elseif options_.lik_init == 4           % Start from the solution of the Riccati equation.
    Pstar = kalman_steady_state(transpose(T),R*Q*transpose(R),transpose(build_selection_matrix(mf,np,vobs)),H);
    Pinf  = [];
    if kalman_algo~=2
        kalman_algo = 1;
    end
elseif options_.lik_init == 5            % Old diffuse Kalman filter only for the non stationary variables
    [eigenvect, eigenv] = eig(T);
    eigenv = diag(eigenv);
    nstable = length(find(abs(abs(eigenv)-1) > 1e-7));
    unstable = find(abs(abs(eigenv)-1) < 1e-7);
    V = eigenvect(:,unstable);
    indx_unstable = find(sum(abs(V),2)>1e-5);
    stable = find(sum(abs(V),2)<1e-5);
    nunit = length(eigenv) - nstable;
    Pstar = options_.Harvey_scale_factor*eye(np);
    if kalman_algo ~= 2
        kalman_algo = 1;
    end
    R_tmp = R(stable, :);
    T_tmp = T(stable,stable);
    Pstar_tmp=lyapunov_solver(T_tmp,R_tmp,Q,DynareOptions);
    Pstar(stable, stable) = Pstar_tmp;
    Pinf  = [];
end
kalman_tol = options_.kalman_tol;
diffuse_kalman_tol = options_.diffuse_kalman_tol;
riccati_tol = options_.riccati_tol;
data1 = Y-trend;
% -----------------------------------------------------------------------------
%  4. Kalman smoother
% -----------------------------------------------------------------------------

if ~missing_value
    for i=1:smpl
        data_index{i}=(1:vobs)';
    end
end

ST = T;
R1 = R;

if options_.heteroskedastic_filter
    Q=get_Qvec_heteroskedastic_filter(Q,smpl,M_);
end

if options_.occbin.smoother.status
    if kalman_algo == 1
        kalman_algo = 2;
    end
    if kalman_algo == 3
        kalman_algo = 4;
    end
end

if kalman_algo == 1 || kalman_algo == 3
    a_initial     = zeros(np,1);
    a_initial=set_Kalman_smoother_starting_values(a_initial,M_,oo_,options_);
    a_initial=T*a_initial; %set state prediction for first Kalman step;
    [alphahat,epsilonhat,etahat,ahat,P,aK,PK,decomp,state_uncertainty, aahat, eehat, d] = missing_DiffuseKalmanSmootherH1_Z(a_initial,ST, ...
        Z,R1,Q,H,Pinf,Pstar, ...
        data1,vobs,np,smpl,data_index, ...
        options_.nk,kalman_tol,diffuse_kalman_tol,options_.filter_decomposition,options_.smoothed_state_uncertainty,options_.filter_covariance,options_.smoother_redux);
    if isinf(alphahat)
        if kalman_algo == 1
            fprintf('\nDsgeSmoother: Switching to univariate filter. This may be a sign of stochastic singularity.\n')
            kalman_algo = 2;
        elseif kalman_algo == 3
            fprintf('\nDsgeSmoother: Switching to univariate filter. This is usually due to co-integration in diffuse filter,\n')
            fprintf(' otherwise it may be a sign of stochastic singularity.\n')
            kalman_algo = 4;
        else
            error('This case shouldn''t happen')
        end
    end
end

if kalman_algo == 2 || kalman_algo == 4
    if ~all(all(abs(H-diag(diag(H)))<1e-14))% ie, the covariance matrix is diagonal...
        if ~expanded_state_vector_for_univariate_filter
            %Augment state vector (follows Section 6.4.3 of DK (2012))
            expanded_state_vector_for_univariate_filter=1;
            Z   = [Z, eye(vobs)];
            ST  = blkdiag(ST,zeros(vobs));
            np  = size(ST,1);
            if options_.heteroskedastic_filter
                Qvec=Q;
                Q=NaN(size(Qvec,1)+size(H,1),size(Qvec,1)+size(H,1),smpl+1);
                for kv=1:size(Qvec,3)
                    Q(:,:,kv) = blkdiag(Qvec(:,:,kv),H);
                end
            else
                Q   = blkdiag(Q,H);
            end
            R1  = blkdiag(R,eye(vobs));
            if kalman_algo == 4
                %recompute Schur state space transformation with
                %expanded state space
                [Pstar,Pinf] = compute_Pinf_Pstar(mf,ST,R1,Q,options_.qz_criterium);
            else
                Pstar = blkdiag(Pstar,H);
                if ~isempty(Pinf)
                    Pinf  = blkdiag(Pinf,zeros(vobs));
                end
            end
            %now reset H to 0
            H   = zeros(vobs,vobs);
        else
            %do nothing, state vector was already expanded
        end
    end
    
    a_initial     = zeros(np,1);
    a_initial=set_Kalman_smoother_starting_values(a_initial,M_,oo_,options_);
    a_initial=ST*a_initial; %set state prediction for first Kalman step;
    [alphahat,epsilonhat,etahat,ahat,P,aK,PK,decomp,state_uncertainty, aahat, eehat, d, regimes_,TT,RR,CC,TTx,RRx,CCx] = missing_DiffuseKalmanSmootherH3_Z(a_initial,ST, ...
        Z,R1,Q,diag(H), ...
        Pinf,Pstar,data1,vobs,np,smpl,data_index, ...
        options_.nk,kalman_tol,diffuse_kalman_tol, ...
        options_.filter_decomposition,options_.smoothed_state_uncertainty,options_.filter_covariance,options_.smoother_redux,occbin_info);
    if options_.occbin.smoother.status
        oo_.occbin.smoother.regime_history = regimes_;
    end
end

if expanded_state_vector_for_univariate_filter && (kalman_algo == 2 || kalman_algo == 4)
    % extracting measurement errors
    % removing observed variables from the state vector
    k = (1:np-vobs);
    alphahat = alphahat(k,:);
    ahat = ahat(k,:);
    aK = aK(:,k,:,:);
    epsilonhat=etahat(end-vobs+1:end,:);
    etahat=etahat(1:end-vobs,:);
    if ~isempty(PK)
        PK = PK(:,k,k,:);
    end
    if ~isempty(decomp)
        decomp = decomp(:,k,:,:);
    end
    if ~isempty(P)
        P = P(k,k,:);
    end
    if ~isempty(state_uncertainty)
        state_uncertainty = state_uncertainty(k,k,:);
    end
end

if ~options_.smoother_redux
    %reset old setting of restricted var_list
    oo_.dr.restrict_var_list = oldoo.restrict_var_list;
    oo_.dr.restrict_columns = oldoo.restrict_columns;
else
    if options_.block == 0
        ic = [ M_.nstatic+(1:M_.nspred) M_.endo_nbr+(1:size(oo_.dr.ghx,2)-M_.nspred) ]';
    else
        ic = oo_.dr.restrict_columns;
    end
    
    if options_.occbin.smoother.status
        % reconstruct occbin smoother
        if length_varargin>0
            % sequence of regimes is provided in input
            isoccbin=1;
        else
            isoccbin=0;
        end
        if length_varargin>1
            TT=varargin{2};
            RR=varargin{3};
            CC=varargin{4};
            if size(TT,3)<(smpl+1)
                TT=repmat(T,1,1,smpl+1);
                RR=repmat(R,1,1,smpl+1);
                CC=repmat(zeros(mm,1),1,smpl+1);
            end
        end
        if isoccbin==0
            [A,B] = kalman_transition_matrix(oo_.dr,(1:M_.endo_nbr)',ic,M_.exo_nbr);
        else
            opts_simul = options_.occbin.simul;
        end
        % reconstruct smoothed variables
        aaa=zeros(M_.endo_nbr,gend);
        aaa(oo_.dr.restrict_var_list,:)=alphahat;
        iTx = zeros(size(TTx));
        for k=1:gend
            if isoccbin
                A = TT(:,:,k);
                B = RR(:,:,k);
                C = CC(:,k);
            else
                C=0;
            end
            iT = pinv(TTx(:,:,k));
            % store pinv
            iTx(:,:,k) = iT;
            Tstar = A(~ismember(1:M_.endo_nbr,oo_.dr.restrict_var_list),oo_.dr.restrict_var_list);
            Rstar = B(~ismember(1:M_.endo_nbr,oo_.dr.restrict_var_list),:);
            Cstar = C(~ismember(1:M_.endo_nbr,oo_.dr.restrict_var_list));
            AS = Tstar*iT;
            BS = Rstar-AS*RRx(:,:,k);
            CS = Cstar-AS*CCx(:,k);
            static_var_list = ~ismember(1:M_.endo_nbr,oo_.dr.restrict_var_list);
            ilagged = any(abs(AS*TTx(:,:,k)-Tstar)'>1.e-12);
            static_var_list0 = static_var_list;
            static_var_list0(static_var_list) = ilagged;
            static_var_list(static_var_list) = ~ilagged;
            aaa(static_var_list,k) = AS(~ilagged,:)*alphahat(:,k)+BS(~ilagged,:)*etahat(:,k)+CS(~ilagged);
            if any(ilagged) && k>1
                aaa(static_var_list0,k) = Tstar(ilagged,:)*alphahat(:,k-1)+Rstar(ilagged,:)*etahat(:,k)+Cstar(ilagged);
            end
            
        end
        alphahat=aaa;

        % reconstruct updated variables
        bbb=zeros(M_.endo_nbr,gend);
        bbb(oo_.dr.restrict_var_list,:)=ahat; % this is t|t
        for k=1:gend
            if isoccbin
                A = TT(:,:,k);
                B = RR(:,:,k);
                C = CC(:,k);
                iT = iTx(:,:,k);
                Tstar = A(~ismember(1:M_.endo_nbr,oo_.dr.restrict_var_list),oo_.dr.restrict_var_list);
                Rstar = B(~ismember(1:M_.endo_nbr,oo_.dr.restrict_var_list),:);
                Cstar = C(~ismember(1:M_.endo_nbr,oo_.dr.restrict_var_list));
                AS = Tstar*iT;
                BS = Rstar-AS*RRx(:,:,k);
                CS = Cstar-AS*CCx(:,k);
                static_var_list = ~ismember(1:M_.endo_nbr,oo_.dr.restrict_var_list);
                ilagged = any(abs(AS*TTx(:,:,k)-Tstar)'>1.e-12);
                static_var_list0 = static_var_list;
                static_var_list0(static_var_list) = ilagged;
                static_var_list(static_var_list) = ~ilagged;
                bbb(static_var_list,k) = AS(~ilagged,:)*ahat(:,k)+BS(~ilagged,:)*eehat(:,k)+CS(~ilagged);
                if any(ilagged) && k>d+1
                    bbb(static_var_list0,k) = Tstar(ilagged,:)*aahat(:,k-1)+Rstar(ilagged,:)*eehat(:,k)+Cstar(ilagged);
                end
            elseif k>d+1
                opts_simul.curb_retrench = options_.occbin.smoother.curb_retrench;
                opts_simul.waitbar = options_.occbin.smoother.waitbar;
                opts_simul.maxit = options_.occbin.smoother.maxit;
                opts_simul.periods = options_.occbin.smoother.periods;
                opts_simul.check_ahead_periods = options_.occbin.smoother.check_ahead_periods;
                opts_simul.full_output = options_.occbin.smoother.full_output;
                opts_simul.piecewise_only = options_.occbin.smoother.piecewise_only;
                opts_simul.SHOCKS = zeros(options_.nk,M_.exo_nbr);
                opts_simul.SHOCKS(1,:) =  eehat(:,k);
                tmp=zeros(M_.endo_nbr,1);
                tmp(oo_.dr.restrict_var_list,1)=aahat(:,k-1);
                opts_simul.endo_init = tmp(oo_.dr.inv_order_var,1);
                opts_simul.init_regime = []; %regimes_(k);
                opts_simul.waitbar=0;
                options_.occbin.simul=opts_simul;
                [~, out] = occbin.solver(M_,oo_,options_);
                % regime in out should be identical to regimes_(k-2) moved one
                % period ahead (so if regimestart was [1 5] it should be [1 4]
                % in out
                %         end
                bbb(oo_.dr.inv_order_var,k) = out.piecewise(1,:) - out.ys';
            end
        end
        ahat0=ahat;
        ahat=bbb;
        if ~isempty(P)
            PP=zeros(M_.endo_nbr,M_.endo_nbr,gend+1);
            PP(oo_.dr.restrict_var_list,oo_.dr.restrict_var_list,:)=P;
            P=PP;
            clear PP
        end
        
        if ~isempty(state_uncertainty)
            mm=size(T,1);
            sstate_uncertainty=zeros(M_.endo_nbr,M_.endo_nbr,gend);
            sstate_uncertainty(oo_.dr.restrict_var_list,oo_.dr.restrict_var_list,:)=state_uncertainty(1:mm,1:mm,:);
            state_uncertainty=sstate_uncertainty;
            clear sstate_uncertainty
        end
        
        aaa = zeros(options_.nk,M_.endo_nbr,gend+options_.nk);
        aaa(:,oo_.dr.restrict_var_list,:)=aK;
        
        for k=2:gend+1
            opts_simul.curb_retrench = options_.occbin.smoother.curb_retrench;
            opts_simul.waitbar = options_.occbin.smoother.waitbar;
            opts_simul.maxit = options_.occbin.smoother.maxit;
            opts_simul.periods = options_.occbin.smoother.periods;
            opts_simul.check_ahead_periods = options_.occbin.smoother.check_ahead_periods;
            opts_simul.full_output = options_.occbin.smoother.full_output;
            opts_simul.piecewise_only = options_.occbin.smoother.piecewise_only;
            opts_simul.SHOCKS = zeros(options_.nk,M_.exo_nbr);
            tmp=zeros(M_.endo_nbr,1);
            tmp(oo_.dr.restrict_var_list,1)=ahat0(:,k-1);
            opts_simul.endo_init = tmp(oo_.dr.inv_order_var,1);
            opts_simul.init_regime = []; %regimes_(k);
            opts_simul.waitbar=0;
            options_.occbin.simul=opts_simul;
            [~, out] = occbin.solver(M_,oo_,options_);
            % regime in out should be identical to regimes_(k-2) moved one
            % period ahead (so if regimestart was [1 5] it should be [1 4]
            % in out
            %         end
            for jnk=1:options_.nk
                aaa(jnk,oo_.dr.inv_order_var,k+jnk-1) = out.piecewise(jnk,:) - out.ys';
            end
        end
        aK=aaa;
        
        if ~isempty(PK)
            PP = zeros(options_.nk,M_.endo_nbr,M_.endo_nbr,gend+options_.nk);
            PP(:,oo_.dr.restrict_var_list,oo_.dr.restrict_var_list,:) = PK;
            PK=PP;
            clear PP
        end
    else
        % reconstruct smoother
        [A,B] = kalman_transition_matrix(oo_.dr,(1:M_.endo_nbr)',ic,M_.exo_nbr);
        iT = pinv(T);
        Tstar = A(~ismember(1:M_.endo_nbr,oo_.dr.restrict_var_list),oo_.dr.restrict_var_list);
        Rstar = B(~ismember(1:M_.endo_nbr,oo_.dr.restrict_var_list),:);
        C = Tstar*iT;
        D = Rstar-C*R;
        static_var_list = ~ismember(1:M_.endo_nbr,oo_.dr.restrict_var_list);
        ilagged = any(abs(C*T-Tstar)'>1.e-12);
        static_var_list0 = static_var_list;
        static_var_list0(static_var_list) = ilagged;
        static_var_list(static_var_list) = ~ilagged;
        % reconstruct smoothed variables
        aaa=zeros(M_.endo_nbr,gend);
        aaa(oo_.dr.restrict_var_list,:)=alphahat;
        for k=1:gend
            aaa(static_var_list,k) = C(~ilagged,:)*alphahat(:,k)+D(~ilagged,:)*etahat(:,k);
        end
        if any(ilagged)
            for k=2:gend
                aaa(static_var_list0,k) = Tstar(ilagged,:)*alphahat(:,k-1)+Rstar(ilagged,:)*etahat(:,k);
            end
        end
        alphahat=aaa;
        
        % reconstruct updated variables
        aaa=zeros(M_.endo_nbr,gend);
        aaa(oo_.dr.restrict_var_list,:)=ahat;
        for k=1:gend
            aaa(static_var_list,k) = C(~ilagged,:)*ahat(:,k)+D(~ilagged,:)*eehat(:,k);
        end
        if any(ilagged)
            %         bbb=zeros(M_.endo_nbr,gend);
            %         bbb(oo_.dr.restrict_var_list,:)=aahat;
            for k=d+2:gend
                aaa(static_var_list0,k) = Tstar(ilagged,:)*aahat(:,k-1)+Rstar(ilagged,:)*eehat(:,k);
            end
        end
        ahat1=aaa;
        % reconstruct aK
        if isempty(options_.nk)
            options_.nk=1;
        end
        aaa = zeros(options_.nk,M_.endo_nbr,gend+options_.nk);
        aaa(:,oo_.dr.restrict_var_list,:)=aK;
        for k=1:gend
            for jnk=1:options_.nk
                aaa(jnk,static_var_list,k+jnk) = C(~ilagged,:)*dynare_squeeze(aK(jnk,:,k+jnk));
            end
        end
        if any(ilagged)
            for k=1:gend
                aaa(1,static_var_list0,k+1) = Tstar(ilagged,:)*ahat(:,k);
                for jnk=2:options_.nk
                    aaa(jnk,static_var_list0,k+jnk) = Tstar(ilagged,:)*dynare_squeeze(aK(jnk-1,:,k+jnk-1));
                end
            end
        end
        aK=aaa;
        ahat=ahat1;
        
        % reconstruct P
        if ~isempty(P)
            PP=zeros(M_.endo_nbr,M_.endo_nbr,gend+1);
            PP(oo_.dr.restrict_var_list,oo_.dr.restrict_var_list,:)=P;
            if ~options_.heteroskedastic_filter                
                DQD=D(~ilagged,:)*Q*transpose(D(~ilagged,:))+C(~ilagged,:)*R*Q*transpose(D(~ilagged,:))+D(~ilagged,:)*Q*transpose(C(~ilagged,:)*R);
                DQR=D(~ilagged,:)*Q*transpose(R);
            end
            for k=1:gend+1
                if options_.heteroskedastic_filter
                    DQD=D(~ilagged,:)*Q(:,:,k)*transpose(D(~ilagged,:))+C(~ilagged,:)*R*Q(:,:,k)*transpose(D(~ilagged,:))+D(~ilagged,:)*Q(:,:,k)*transpose(C(~ilagged,:)*R);
                    DQR=D(~ilagged,:)*Q(:,:,k)*transpose(R);
                end
                PP(static_var_list,static_var_list,k)=C(~ilagged,:)*P(:,:,k)*C(~ilagged,:)'+DQD;
                PP(static_var_list,oo_.dr.restrict_var_list,k)=C(~ilagged,:)*P(:,:,k)+DQR;
                PP(oo_.dr.restrict_var_list,static_var_list,k)=transpose(PP(static_var_list,oo_.dr.restrict_var_list,k));
            end
            P=PP;
            clear PP
        end
        
        % reconstruct state_uncertainty
        if ~isempty(state_uncertainty)
            mm=size(T,1);
            ss=length(find(static_var_list));
            sstate_uncertainty=zeros(M_.endo_nbr,M_.endo_nbr,gend);
            sstate_uncertainty(oo_.dr.restrict_var_list,oo_.dr.restrict_var_list,:)=state_uncertainty(1:mm,1:mm,:);
            for k=1:gend
                sstate_uncertainty(static_var_list,static_var_list,k)=[C(~ilagged,:) D(~ilagged,:)]*state_uncertainty(:,:,k)*[C(~ilagged,:) D(~ilagged,:)]';
                tmp = [C(~ilagged,:) D(~ilagged,:)]*state_uncertainty(:,:,k);
                sstate_uncertainty(static_var_list,oo_.dr.restrict_var_list,k)=tmp(1:ss,1:mm);
                sstate_uncertainty(oo_.dr.restrict_var_list,static_var_list,k)=transpose(sstate_uncertainty(static_var_list,oo_.dr.restrict_var_list,k));
            end
            state_uncertainty=sstate_uncertainty;
            clear sstate_uncertainty
        end
        
        % reconstruct PK
        if ~isempty(PK)
            PP = zeros(options_.nk,M_.endo_nbr,M_.endo_nbr,gend+options_.nk);
            PP(:,oo_.dr.restrict_var_list,oo_.dr.restrict_var_list,:) = PK;
            if ~options_.heteroskedastic_filter
                DQD=D(~ilagged,:)*Q*transpose(D(~ilagged,:))+C(~ilagged,:)*R*Q*transpose(D(~ilagged,:))+D(~ilagged,:)*Q*transpose(C(~ilagged,:)*R);
                DQR=D(~ilagged,:)*Q*transpose(R);
                for f=1:options_.nk
                    for k=1:gend
                        PP(f,static_var_list,static_var_list,k+f)=C(~ilagged,:)*squeeze(PK(f,:,:,k+f))*C(~ilagged,:)'+DQD;
                        PP(f,static_var_list,oo_.dr.restrict_var_list,k+f)=C(~ilagged,:)*squeeze(PK(f,:,:,k+f))+DQR;
                        PP(f,oo_.dr.restrict_var_list,static_var_list,k+f)=transpose(squeeze(PP(f,static_var_list,oo_.dr.restrict_var_list,k+f)));
                    end
                end
            end
            PK=PP;
            clear PP
        end
    end
    
    bayestopt_.mf = bayestopt_.smoother_var_list(bayestopt_.smoother_mf);
end

function a=set_Kalman_smoother_starting_values(a,M_,oo_,options_)
% function a=set_Kalman_smoother_starting_values(a,M_,oo_,options_)
% Sets initial states guess for Kalman filter/smoother based on M_.filter_initial_state
%
% INPUTS
%   o a             [double]   (p*1) vector of states
%   o M_            [structure] decribing the model
%   o oo_           [structure] storing the results
%   o options_      [structure] describing the options
%
% OUTPUTS
%   o a             [double]    (p*1) vector of set initial states

if isfield(M_,'filter_initial_state') && ~isempty(M_.filter_initial_state)
    state_indices=oo_.dr.order_var(oo_.dr.restrict_columns);
    for ii=1:size(state_indices,1)
        if ~isempty(M_.filter_initial_state{state_indices(ii),1})
            if options_.loglinear && ~options_.logged_steady_state
                a(oo_.dr.restrict_columns(ii)) = log(eval(M_.filter_initial_state{state_indices(ii),2})) - log(oo_.dr.ys(state_indices(ii)));
            elseif ~options_.loglinear && ~options_.logged_steady_state
                a(oo_.dr.restrict_columns(ii)) = eval(M_.filter_initial_state{state_indices(ii),2}) - oo_.dr.ys(state_indices(ii));
            else
                error(['The steady state is logged. This should not happen. Please contact the developers'])
            end
        end
    end
end

