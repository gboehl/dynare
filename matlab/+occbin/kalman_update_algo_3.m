function [a, a1, P, P1, v, Fi, Ki, T, R, C, regimes_, error_flag, M_, lik, alphahat, etahat, TT, RR, CC] = kalman_update_algo_3(a,a1,P,P1,data_index,Z,v,Fi,Ki,Y,H,QQQ,T0,R0,TT,RR,CC,regimes0,M_,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,options_,occbin_options,kalman_tol,nk)
% function [a, a1, P, P1, v, Fi, Ki, T, R, C, regimes_, error_flag, M_, alphahat, etahat, TT, RR, CC] = kalman_update_algo_3(a,a1,P,P1,data_index,Z,v,Fi,Ki,Y,H,QQQ,T0,R0,TT,RR,CC,regimes0,M_,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,options_,occbin_options,kalman_tol,nk)
%
% INPUTS
% - a               [N by 1]                t-1's state estimate
% - a1              [N by N by 2]           state predictions made at t-1:t
% - P               [N by N]                t-1's covariance of states
% - P1              [N by N by 2]           one-step ahead forecast error variance at t-1:t
% - data_index:     [cell]                  1*2 cell of column vectors of indices.
% - Z               [N_obs ny N]            Selector matrix
% - v               [N_obs by 2]            prediction error on observables at t-1:t
% - Fi              [N_obs by 1]            F_i matrix
% - Ki              [N by N_obs]            Kalman gain matrix
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
% - options_        [structure]             Matlab's structure describing the current options (options_).
% - dr                   [structure]        Reduced form model.
% - endo_steady_state    [vector]           steady state value for endogenous variables
% - exo_steady_state     [vector]           steady state value for exogenous variables
% - exo_det_steady_state [vector]           steady state value for exogenous deterministic variables
% - occbin_options_ [structure]             Matlab's structure describing the Occbin options.
% - kalman_tol      [double]                tolerance for reciprocal condition number
% - nk              [double]                number of forecasting periods
% 
% Outputs
% - a               [N by 2]                t-1's state estimate
% - a1              [N by N by 2]           state predictions made at t-1:t
% - P               [N by N by 2]           t-1's covariance of states
% - P1              [N by N by 2]           one-step ahead forecast error variance at t-1:t
% - v               [N_obs by 2]            prediction error on observables at t-1:t
% - Fi              [N_obs by 2]            F_i matrix
% - Ki              [N by N_obs by 2]       Kalman gain matrix
% - TT              [N by N by 2]           state transition matrix at t-1:t
% - RR              [N by N_exo by 2]       shock impact matrix at t-1:t
% - CC              [N by 2]                state space constant state transition matrix at t-1:t
% - regimes_        [structure]             regime info at t-1:t
% - error_flag      [structure]             error flag
% - M_              [structure]             Matlab's structure describing the model (M_).
% - alphahat:                               smoothed variables (a_{t|T})
% - etahat:                                 smoothed shocks
% - TT              [N by N by 2]           state transition matrix at t-1:t
% - RR              [N by N_exo by 2]       shock impact matrix at t-1:t
% - CC              [N by 2]                state space constant state transition matrix at t-1:t
%
% Notes: The algorithm and implementation is based on Massimo Giovannini,
% Philipp Pfeiffer, Marco Ratto (2021), Efficient and robust inference of models with occasionally binding
% constraints, Working Papers 2021-03, Joint Research Centre, European Commission 


% Copyright © 2021-2024 Dynare Team
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

warning off

options_.noprint = true;
T=[];
R=[];
C=[];
regimes_ = struct();

if isempty(nk)
    nk=1;
end
nk=max(nk,1);

opts_simul = occbin_options.opts_simul;
options_.occbin.simul=opts_simul;
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
myrestrict=[];
if options_.smoother_redux
    opts_simul.restrict_state_space =1;
    myrestrict='restrict';
end

mm=size(a,1);
if ~isempty(fieldnames(regimes0))
    if options_.occbin.filter.guess_regime
        [~,~,~,~,~,~, TTx, RRx, CCx] ...
            = occbin.dynare_resolve(M_,options_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state, regimes0(1),myrestrict,T0,R0);
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
%         regimes0=[regimes0 base_regime base_regime];
        TT(:,:,2) = TTx(:,:,end);
        RR(:,:,2) = RRx(:,:,end);
        CC(:,2) = CCx(:,end);
    end
    [a, a1, P, P1, v, Fi, Ki, alphahat, etahat, lik] = occbin_kalman_update(a,a1,P,P1,data_index,Z,v,Y,H,QQQ,TT,RR,CC,Ki,Fi,mm,kalman_tol);
else
    [~,~,~,~,~,~, TTx, RRx, CCx] ...
        = occbin.dynare_resolve(M_,options_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state, base_regime,myrestrict,T0,R0);
    if isempty(fieldnames(regimes0))
        regimes0 = regimes_;
    else
        regimes0(1)=base_regime;
    end
    TT(:,:,2) = TTx(:,:,end);
    RR(:,:,2) = RRx(:,:,end);
    CC(:,2) = CCx(:,end);
    [a, a1, P, P1, v, Fi, Ki, alphahat, etahat, lik] = occbin_kalman_update(a,a1,P,P1,data_index,Z,v,Y,H,QQQ,TT,RR,CC,Ki,Fi,mm,kalman_tol);
    regimes0(1)=base_regime;
end


%% run here the occbin simul
opts_simul.SHOCKS = zeros(max(3,1+nk),M_.exo_nbr);
opts_simul.SHOCKS(1,:) = etahat(:,2)';
opts_simul.exo_pos=1:M_.exo_nbr;

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
    return;
end

regimes_ = out.regime_history;
if M_.occbin.constraint_nbr==1
    myregime = [regimes_.regime];
else
    myregime = [regimes_.regime1 regimes_.regime2];
end
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
        [a, a1, P, P1, v, Fi, Ki, alphahat, etahat, lik] = occbin_kalman_update(a,a1,P,P1,data_index,Z,v,Y,H,QQQ,TT,RR,CC,Ki,Fi,mm,kalman_tol);
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
            error_flag = out.error_flag;
            etahat=etahat(:,2);
            lik=inf;
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
                        return;
                    else
                        regimes_ = out.regime_history;
                        TT(:,:,2)=ss.T(my_order_var,my_order_var,1);
                        RR(:,:,2)=ss.R(my_order_var,:,1);
                        CC(:,2)=ss.C(my_order_var,1);
                        [a, a1, P, P1, v, Fi, Ki, alphahat, etahat, lik] = occbin_kalman_update(a,a1,P,P1,data_index,Z,v,Y,H,QQQ,TT,RR,CC,Ki,Fi,mm,kalman_tol);
                    end
                    else
                        error_flag = 330;
                        etahat=etahat(:,2);
                        lik=inf;
                        return;
                    end
            end
        end
    end
end

error_flag = out.error_flag;
if ~error_flag && niter>options_.occbin.likelihood.max_number_of_iterations && ~isequal(regimes_(1),regimes0(1)) %fixed point algorithm did not converge
  error_flag = 1;
end

if ~error_flag
    a = out.piecewise(1:nk+1,my_order_var)' - repmat(out.ys(my_order_var),1,nk+1);
    regimes_=regimes_(1:3);
end
T = ss.T(my_order_var,my_order_var,:);
R = ss.R(my_order_var,:,:);
C = ss.C(my_order_var,:);
TT = ss.T(dr.order_var,dr.order_var,1);
RR = ss.R(dr.order_var,:,1);
CC = ss.C(dr.order_var,1);
QQ = R(:,:,2)*QQQ(:,:,3)*transpose(R(:,:,2));
P(:,:,1) = P(:,:,2);
for j=1:nk
    P(:,:,j+1) = T(:,:,j+1)*P(:,:,j)*transpose(T(:,:,j+1))+QQ;
end
% P = cat(3,P(:,:,2),P2);


end

function [a, a1, P, P1, v, Fi, Ki, alphahat, etahat, lik] = occbin_kalman_update(a,a1,P,P1,data_index,Z,v,Y,H,QQQ,TT,RR,CC,Ki,Fi,mm,kalman_tol)
% [a, a1, P, P1, v, Fi, Ki, alphahat, etahat] = occbin_kalman_update(a,a1,P,P1,data_index,Z,v,Y,H,QQQ,TT,RR,CC,Ki,Fi,mm,kalman_tol)
% - a
% - a1
% - P
% - P1
% - v
% - Fi
% - Ki


t=2;

%% forward pass
% given updated variables and covariance in t=1, we make the step to t=2
T = TT(:,:,t);
R = RR(:,:,t);
C = CC(:,t);
Q=QQQ(:,:,t);
QQ = R*Q*transpose(R);
a1(:,t) = T*a(:,t-1)+C;
a(:,t) = a1(:,t);
P1(:,:,t) = T*P(:,:,t-1)*T' + QQ;                                        %transition according to (6.14) in DK (2012)
P(:,:,t) = P1(:,:,t);
di = data_index{t}';
% store info for lik
if not(isempty(di))
    vv      = Y(di,t) - Z(di,:)*a(:,t);
    F = Z(di,:)*P(:,:,t)*Z(di,:)' + diag(H(di));
    sig=sqrt(diag(F));
    lik=inf;
end

if isempty(di)
    Fi(:,t) = 0;
    Ki(:,:,t) = 0;
    lik =0;
end
for i=di
    Zi = Z(i,:);
    v(i,t)  = Y(i,t) - Zi*a(:,t);                                       % nu_{t,i} in 6.13 in DK (2012)
    Fi(i,t) = Zi*P(:,:,t)*Zi' + H(i);                                   % F_{t,i} in 6.13 in DK (2012), relies on H being diagonal
    Ki(:,i,t) = P(:,:,t)*Zi';                                           % K_{t,i}*F_(i,t) in 6.13 in DK (2012)
    if Fi(i,t) > kalman_tol
        a(:,t) = a(:,t) + Ki(:,i,t)*v(i,t)/Fi(i,t);                     %filtering according to (6.13) in DK (2012)
        P(:,:,t) = P(:,:,t) - Ki(:,i,t)*Ki(:,i,t)'/Fi(i,t);             %filtering according to (6.13) in DK (2012)
    else
        % do nothing as a_{t,i+1}=a_{t,i} and P_{t,i+1}=P_{t,i}, see
        % p. 157, DK (2012)
    end
end
if not(isempty(di))
    log_dF = log(det(F./(sig*sig')))+2*sum(log(sig));
    iF = inv(F./(sig*sig'))./(sig*sig');
    lik = log_dF + transpose(vv)*iF*vv + length(di)*log(2*pi);
end

%% do backward pass
ri=zeros(mm,1);
t = t+1;
while t > 1
    t = t-1;
    di = flipud(data_index{t})';
    for i = di
        if Fi(i,t) > kalman_tol
            Li = eye(mm)-Ki(:,i,t)*Z(i,:)/Fi(i,t);
            ri = Z(i,:)'/Fi(i,t)*v(i,t)+Li'*ri;                             % DK (2012), 6.15, equation for r_{t,i-1}
        end
    end
    r(:,t) = ri;                                                            % DK (2012), below 6.15, r_{t-1}=r_{t,0}
    alphahat(:,t) = a1(:,t) + P1(:,:,t)*r(:,t);
    Q=QQQ(:,:,t);
    QRt = Q*transpose(RR(:,:,t));
    T = TT(:,:,t);
    etahat(:,t) = QRt*r(:,t);
    ri = T'*ri;                                                             % KD (2003), eq. (23), equation for r_{t-1,p_{t-1}}
end

end
