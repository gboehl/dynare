function [a, a1, P, P1, v, T, R, C, regimes_, error_flag, M_, lik, etahat, alphahat, V] = kalman_update_algo_1(a,a1,P,P1,data_index,Z,v,Y,H,QQQ,T0,R0,TT,RR,CC,regimes0,M_,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,options_,occbin_options)
% function [a, a1, P, P1, v, T, R, C, regimes_, error_flag, M_, lik, etahat] = kalman_update_algo_1(a,a1,P,P1,data_index,Z,v,Y,H,QQQ,T0,R0,TT,RR,CC,regimes0,M_,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,options_,occbin_options)
% INPUTS
% - a               [N by 1]                t-1's state estimate
% - a1              [N by 2]               state predictions made at t-1:t
% - P               [N by N]                t-1's covariance of states
% - P1              [N by N by 2]           one-step ahead forecast error variance at t-1:t
% - data_index:     [cell]                  1*2 cell of column vectors of indices.
% - Z               [N_obs ny N]            Selector matrix
% - v               [N_obs by 2]            prediction error on observables at t-1:t
% - Y:              [N_obs by 2]            observations at t-1:t
% - H               [N_obs by 1]            vector of measurement error
% - QQQ             [N_exo by N_exo by 3]   covariance matrix of shocks at t-1:t+1
% - T0              [N by N]                initial state transition matrix
% - R0              [N by N_exo]            initial shock impact transition matrix
% - TT              [N by N by 2]           state transition matrix at t-1:t
% - RR              [N by N_exo by 2]       shock impact matrix at t-1:t
% - CC              [N by 2]                state space constant state transition matrix at t-1:t
% - regimes0        [structure]             regime info at t-1:t
% - M_              [structure]             Matlab's structure describing the model (M_).
% - dr                   [structure]        Reduced form model.
% - endo_steady_state    [vector]           steady state value for endogenous variables
% - exo_steady_state     [vector]           steady state value for exogenous variables
% - exo_det_steady_state [vector]           steady state value for exogenous deterministic variables
% - options_        [structure]             Matlab's structure describing the current options (options_).
% - occbin_options_ [structure]             Matlab's structure describing the Occbin options.
% - kalman_tol      [double]                tolerance for reciprocal condition number
% 
% Outputs
% - a               [N by 2]                t-1's state estimate
% - a1              [N by 2]                state predictions made at t-1:t
% - P               [N by N by 2]           t-1's covariance of states
% - P1              [N by N by 2]           one-step ahead forecast error variance at t-1:t
% - v               [N_obs by 2]            prediction error on observables at t-1:t
% - T               [N by N by 2]           state transition matrix at t-1:t
% - R               [N by N_exo by 2]       shock impact matrix at t-1:t
% - C               [N by 2]                state space constant state transition matrix at t-1:t
% - regimes_        [structure]             regime info at t-1:t
% - error_flag      [integer]               error code
% - M_              [structure]             Matlab's structure describing the model (M_).
% - lik             [double]                likelihood
% - etahat:                                 smoothed shocks
%
% Notes: The algorithm and implementation is based on Massimo Giovannini,
% Philipp Pfeiffer, Marco Ratto (2021), Efficient and robust inference of models with occasionally binding
% constraints, Working Papers 2021-03, Joint Research Centre, European Commission 

% Copyright © 2021-2023 Dynare Team
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

orig_warning_state = warning;

options_.noprint = true;
R=NaN(size(RR));
C=NaN(size(CC));
T=NaN(size(TT));
lik=Inf;

sto.a=a;
sto.a1=a1;
sto.P=P;
sto.P1=P1;

base_regime = struct();
if M_.occbin.constraint_nbr==1
    base_regime.regime = 0;
    base_regime.regimestart = 1;
else
    base_regime.regime1 = 0;
    base_regime.regimestart1 = 1;
    base_regime.regime2 = 0;
    base_regime.regimestart2 = 1;
end
regimes_ = [base_regime base_regime base_regime];
opts_simul = occbin_options.opts_simul;
options_.occbin.simul=opts_simul;

mm=size(a,1);
%% store info in t=1
t=1;
di = data_index{t};
T = TT(:,:,t);
ZZ = Z(di,:);
di = data_index{t};
F = ZZ*P1(:,:,t)*ZZ' + H(di,di);
Fi(di,di,t)=F;
sig=sqrt(diag(F));  
iF(di,di,t)   = inv(F./(sig*sig'))./(sig*sig');
PZI         = P1(:,:,t)*ZZ'*iF(di,di,t);
% K(:,di,t)    = T*PZI;
% L(:,:,t)    = T-K(:,di,t)*ZZ;
L(:,:,t)    = eye(mm)-PZI*ZZ;

if ~isempty(fieldnames(regimes0))
    if options_.occbin.filter.guess_regime
        [~,~,~,~,~,~, TTx, RRx, CCx] ...
            = occbin.dynare_resolve(M_,options_,dr, endo_steady_state, exo_steady_state, exo_det_steady_state, regimes0(1),'reduced_state_space',T0,R0);
        if M_.occbin.constraint_nbr==1
            bindx = occbin.backward_map_regime(regimes0(1).regime, regimes0(1).regimestart);
            bindx = bindx(2:end);
            [regimes0(2).regime, regimes0(2).regimestart, error_flag]=occbin.map_regime(bindx);
            bindx = bindx(2:end);
            [regimes0(3).regime, regimes0(3).regimestart, error_flag]=occbin.map_regime(bindx);
        else
            bindx1 = occbin.backward_map_regime(regimes0(1).regime1, regimes0(1).regimestart1);
            bindx2 = occbin.backward_map_regime(regimes0(1).regime2, regimes0(1).regimestart2);
            bindx1 = bindx1(2:end);
            bindx2 = bindx2(2:end);
            [regimes0(2).regime1, regimes0(2).regimestart1, error_flag]=occbin.map_regime(bindx1);
            [regimes0(2).regime2, regimes0(2).regimestart2, error_flag]=occbin.map_regime(bindx2);
            bindx1 = bindx1(2:end);
            bindx2 = bindx2(2:end);
            [regimes0(3).regime1, regimes0(3).regimestart1, error_flag]=occbin.map_regime(bindx1);
            [regimes0(3).regime2, regimes0(3).regimestart2, error_flag]=occbin.map_regime(bindx2);
        end
%     regimes0=[regimes0 base_regime base_regime];
        TT(:,:,2) = TTx(:,:,end);
        RR(:,:,2) = RRx(:,:,end);
        CC(:,2) = CCx(:,end);
    end
    [a, a1, P, P1, v, alphahat, etahat, lik, V, error_flag] = occbin_kalman_update0(a,a1,P,P1,data_index,Z,v,Y,H,QQQ,TT,RR,CC,iF,L,mm, options_.rescale_prediction_error_covariance, options_.occbin.likelihood.IF_likelihood, options_.occbin.filter.state_covariance);
else
    [~,~,~,~,~,~, TTx, RRx, CCx] ...
        = occbin.dynare_resolve(M_,options_,dr, endo_steady_state, exo_steady_state, exo_det_steady_state, base_regime,'reduced_state_space',T0,R0);
    if isempty(fieldnames(regimes0))
        regimes0 = regimes_;
    else
        regimes0(1)=base_regime;
    end
    TT(:,:,2) = TTx(:,:,end);
    RR(:,:,2) = RRx(:,:,end);
    CC(:,2) = CCx(:,end);
    [a, a1, P, P1, v, alphahat, etahat, lik, V, error_flag] = occbin_kalman_update0(a,a1,P,P1,data_index,Z,v,Y,H,QQQ,TT,RR,CC,iF,L,mm, options_.rescale_prediction_error_covariance, options_.occbin.likelihood.IF_likelihood, options_.occbin.filter.state_covariance);
end
if error_flag
    etahat=NaN(size(QQQ,1),1);
    warning(orig_warning_state);
    return;
end

%% run here the occbin simul
opts_simul.SHOCKS = zeros(3,M_.exo_nbr);
opts_simul.exo_pos=1:M_.exo_nbr;
opts_simul.SHOCKS(1,:) = etahat(:,2)';
if opts_simul.restrict_state_space
    tmp=zeros(M_.endo_nbr,1);
    tmp(dr.restrict_var_list,1)=alphahat(:,1);
    opts_simul.endo_init = tmp(dr.inv_order_var,1);
    my_order_var = dr.order_var(dr.restrict_var_list);
else
    opts_simul.endo_init = alphahat(dr.inv_order_var,1);
    my_order_var = dr.order_var;
end
options_.occbin.simul=opts_simul;
if options_.occbin.filter.guess_regime
    options_.occbin.simul.init_regime=regimes0;
    [~, out, ss] = occbin.solver(M_,options_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state);
    if out.error_flag
        options_.occbin.simul=opts_simul;
        [~, out, ss] = occbin.solver(M_,options_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state);
    end
else
    [~, out, ss] = occbin.solver(M_,options_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state);
    if out.error_flag
        options_.occbin.simul.init_regime=regimes0;
        [~, out, ss] = occbin.solver(M_,options_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state);
    end
end
if out.error_flag
    error_flag = out.error_flag;
    etahat=etahat(:,2);
    lik=inf;
    warning(orig_warning_state);
    return;
end

regimes_ = out.regime_history;
if M_.occbin.constraint_nbr==1
    myregime = [regimes_.regime];
else
    myregime = [regimes_.regime1 regimes_.regime2];
end
etahat_hist = {etahat};
regime_hist = {regimes0(1)};
if M_.occbin.constraint_nbr==1
    regime_end = regimes0(1).regimestart(end);
end
lik_hist=lik;
niter=1;
is_periodic=0;

if any(myregime) || ~isequal(regimes_(1),regimes0(1))
    while ~isequal(regimes_(1),regimes0(1)) && ~is_periodic && ~out.error_flag && niter<=options_.occbin.likelihood.max_number_of_iterations
        niter=niter+1;
        TT(:,:,2)=ss.T(my_order_var,my_order_var,1);
        RR(:,:,2)=ss.R(my_order_var,:,1);
        CC(:,2)=ss.C(my_order_var,1);
        regime_hist(niter) = {regimes_(1)};
        if M_.occbin.constraint_nbr==1
            regime_end(niter) = regimes_(1).regimestart(end);
        end
        [a, a1, P, P1, v, alphahat, etahat, lik, V] = occbin_kalman_update0(a,a1,P,P1,data_index,Z,v,Y,H,QQQ,TT,RR,CC,iF,L,mm, options_.rescale_prediction_error_covariance, options_.occbin.likelihood.IF_likelihood, options_.occbin.filter.state_covariance);
        etahat_hist(niter) = {etahat};
        lik_hist(niter) = lik;
        opts_simul.SHOCKS(1,:) = etahat(:,2)';
        if opts_simul.restrict_state_space
            tmp=zeros(M_.endo_nbr,1);
            tmp(dr.restrict_var_list,1)=alphahat(:,1);
            opts_simul.endo_init = tmp(dr.inv_order_var,1);
        else
            opts_simul.endo_init = alphahat(dr.inv_order_var,1);
        end
        opts_simul.init_regime=regimes_(1);
        if M_.occbin.constraint_nbr==1
            myregimestart = [regimes_.regimestart];
        else
            myregimestart = [regimes_.regimestart1 regimes_.regimestart2];
        end
        opts_simul.periods = max(opts_simul.periods,max(myregimestart));
        options_.occbin.simul=opts_simul;
        [~, out, ss] = occbin.solver(M_,options_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state);
        if out.error_flag
            options_.occbin.simul.init_regime=[];
            [~, out, ss] = occbin.solver(M_,options_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state);
        end
        if out.error_flag
            error_flag = out.error_flag;
            etahat=etahat(:,2);
            lik=inf;
            warning(orig_warning_state);
            return;
        end
        regimes0=regimes_;
        regimes_ = out.regime_history;
        if niter>1
            for kiter=1:niter-1
                is_periodic(kiter) = isequal(regime_hist{kiter}, regimes_(1));
            end
            is_periodic_iter = find(is_periodic);
            is_periodic =  any(is_periodic);
            if is_periodic
                % re-set to previous regime
                if options_.occbin.filter.periodic_solution
                    % force projection conditional on most likely regime
                    [m, im]=min(lik_hist(is_periodic_iter:end));
                    opts_simul.init_regime=regime_hist{is_periodic_iter+im-1};
                    if M_.occbin.constraint_nbr==1
                        myregimestart = [regimes0.regimestart];
                    else
                        myregimestart = [regimes0.regimestart1 regimes0.regimestart2];
                    end
                    opts_simul.periods = max(opts_simul.periods,max(myregimestart));
                    opts_simul.maxit=1;
                    options_.occbin.simul=opts_simul;
                    [~, out, ss] = occbin.solver(M_,options_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state);
                    if out.error_flag
                        error_flag = out.error_flag;
                        etahat=etahat(:,2);
                        lik=inf;
                        warning(orig_warning_state);
                        return;
                    else
                        regimes_ = out.regime_history;
                        TT(:,:,2)=ss.T(my_order_var,my_order_var,1);
                        RR(:,:,2)=ss.R(my_order_var,:,1);
                        CC(:,2)=ss.C(my_order_var,1);
                        [a, a1, P, P1, v, alphahat, etahat, lik, V] = occbin_kalman_update0(a,a1,P,P1,data_index,Z,v,Y,H,QQQ,TT,RR,CC,iF,L,mm, options_.rescale_prediction_error_covariance, options_.occbin.likelihood.IF_likelihood, options_.occbin.filter.state_covariance);
                    end
                else
                    error_flag = 330;
                    etahat=etahat(:,2);
                    lik=inf;
                    warning(orig_warning_state);
                    return;
                end
            end
        end
    end
end

error_flag = out.error_flag;
if ~error_flag && niter>options_.occbin.likelihood.max_number_of_iterations && ~isequal(regimes_(1),regimes0(1))
    error_flag = 1;
end

if ~error_flag
    a = out.piecewise(1:2,my_order_var)' - repmat(out.ys(my_order_var),1,2);
    regimes_=regimes_(1:3);
end
T = ss.T(my_order_var,my_order_var,1:2);
R = ss.R(my_order_var,:,1:2);
C = ss.C(my_order_var,1:2);
QQ = R(:,:,2)*QQQ(:,:,3)*transpose(R(:,:,2));
P(:,:,1) = P(:,:,2);
P(:,:,2) = T(:,:,2)*P(:,:,1)*transpose(T(:,:,2))+QQ;
if nargout>=13
    etahat=etahat(:,2);
end

%reset warning state
warning(orig_warning_state);
end

function [a, a1, P, P1, v, alphahat, etahat, lik, V, error_flag] = occbin_kalman_update0(a,a1,P,P1,data_index,Z,v,Y,H,QQQ,TT,RR,CC,iF,L,mm, rescale_prediction_error_covariance, IF_likelihood, state_uncertainty_flag)
alphahat=NaN(size(a));
etahat=NaN(size(QQQ,1),2);
lik=Inf; 
error_flag=0;

if state_uncertainty_flag
    V               = zeros(mm,mm,2);
    N               = zeros(mm,mm,3);
else
    V=[];
end

orig_warning_state = warning;
if nargin<18
    IF_likelihood=0;
end
t=2;
lik=0;
%% forward pass
% given updated variables and covarnace in t=1, we make the step to t=2
T = TT(:,:,t);
R = RR(:,:,t);
C = CC(:,t);
Q=QQQ(:,:,t);
QQ = R*Q*transpose(R);
a1(:,t) = T*a(:,t-1)+C;
a(:,t) = a1(:,t);
P1(:,:,t) = T*P(:,:,t-1)*T' + QQ;                                        %transition according to (6.14) in DK (2012)
P(:,:,t) = P1(:,:,t);

di = data_index{t};
if isempty(di)
    a(:,t)     = a1(:,t);
    L(:,:,t)        = eye(mm);
    P1(:,:,t+1)      = T*P(:,:,t)*T' + QQ;                               %p. 111, DK(2012)
else
    ZZ = Z(di,:);
    v(di,t)      = Y(di,t) - ZZ*a(:,t);
    F = ZZ*P(:,:,t)*ZZ' + H(di,di);
    sig=sqrt(diag(F));
    if any(any(isnan(F)))
        error_flag=1;
        warning(orig_warning_state);
        return;
    end
    if rank(F)<size(F,1) 
        % here we trap cases when some OBC regime triggers singularity 
        % e.g. no shock to interest rate at ZLB
        di=di(find(sig));
        ZZ = Z(di,:);
        v(di,t)      = Y(di,t) - ZZ*a(:,t);
        F = ZZ*P(:,:,t)*ZZ' + H(di,di);
        sig=sqrt(diag(F));
        data_index{t}=di;
    end

    if rescale_prediction_error_covariance
        if IF_likelihood == 0
            log_dF = log(det(F./(sig*sig')))+2*sum(log(sig));
        end
        iF(di,di,t) = inv(F./(sig*sig'))./(sig*sig');
    else
        if IF_likelihood == 0
            log_dF = log(det(F));
        end
        iF(di,di,t) = inv(F);
    end
    if IF_likelihood == 0
        lik = log_dF + transpose(v(di,t))*iF(di,di,t)*v(di,t) + length(di)*log(2*pi);
    end
    
    PZI         = P(:,:,t)*ZZ'*iF(di,di,t);
    a(:,t) = a(:,t) + PZI*v(di,t);
    K    = PZI*ZZ;
    L(:,:,t)    = (eye(mm) - K);
    P(:,:,t)  = P(:,:,t) - K*P(:,:,t);
    
end

%% do backward pass
r=zeros(mm,3);
t = t+1;
while t > 1
    t = t-1;
    di = data_index{t};
    if isempty(di)
        % in this case, L is simply T due to Z=0, so that DK (2012), eq. 4.93 obtains
        r(:,t) = L(:,:,t)'*r(:,t+1);                                        %compute r_{t-1}, DK (2012), eq. 4.38 with Z=0
        if state_uncertainty_flag
            N(:,:,t)=L(:,:,t)'*N(:,:,t+1)*L(:,:,t);                         %compute N_{t-1}, DK (2012), eq. 4.42 with Z=0
        end
    else
        ZZ = Z(di,:);
        r(:,t) = ZZ'*iF(di,di,t)*v(di,t) + L(:,:,t)'*r(:,t+1);              %compute r_{t-1}, DK (2012), eq. 4.38
        if state_uncertainty_flag
            N(:,:,t)=ZZ'*iF(di,di,t)*ZZ+L(:,:,t)'*N(:,:,t+1)*L(:,:,t);      %compute N_{t-1}, DK (2012), eq. 4.42
        end
    end
    Q=QQQ(:,:,t);
    QRt = Q*transpose(RR(:,:,t));
    T = TT(:,:,t);
    alphahat(:,t)       = a1(:,t) + P1(:,:,t)*r(:,t);                         %DK (2012), eq. 4.35
    etahat(:,t) = QRt*r(:,t);                                               %DK (2012), eq. 4.63
    r(:,t) = T'*r(:,t);                                                             % KD (2003), eq. (23), equation for r_{t-1,p_{t-1}}
    if state_uncertainty_flag
        V(:,:,t)    = P1(:,:,t)-P1(:,:,t)*N(:,:,t)*P1(:,:,t);                      %DK (2012), eq. 4.43
        N(:,:,t)    = T'*N(:,:,t)*T;                                         %DK (2012), eq. 4.43
    end

    if IF_likelihood && t==2 && not(isempty(di))
        ishocks = any(ZZ*RR(:,:,t));
        ishocks(find(etahat(ishocks,t)==0))=false;
        Gmat1  = ZZ*RR(:,ishocks,t);
        if size(Gmat1,1) == size(Gmat1,2)
            log_det_jacobian = log(det(Q(ishocks,ishocks))) + 2*log(abs(det(Gmat1)));
            trace_term = etahat(ishocks,t)'*(Q(ishocks,ishocks)\etahat(ishocks,t));
            
            lik = log_det_jacobian + trace_term  + length(di)*log(2*pi);
        else
            lik = inf;
        end
        
    end
end

%reset warning state
warning(orig_warning_state);
end
