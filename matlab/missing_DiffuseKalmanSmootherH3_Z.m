function [alphahat,epsilonhat,etahat,a,P1,aK,PK,decomp,V, aalphahat,eetahat,d,varargout] = missing_DiffuseKalmanSmootherH3_Z(a_initial,T,Z,R,Q,H,Pinf1,Pstar1,Y,pp,mm,smpl,data_index,nk,kalman_tol,diffuse_kalman_tol,decomp_flag,state_uncertainty_flag, filter_covariance_flag, smoother_redux, occbin_)
% function [alphahat,epsilonhat,etahat,a,P1,aK,PK,decomp,V, aalphahat,eetahat,d] = missing_DiffuseKalmanSmootherH3_Z(a_initial,T,Z,R,Q,H,Pinf1,Pstar1,Y,pp,mm,smpl,data_index,nk,kalman_tol,diffuse_kalman_tol,decomp_flag,state_uncertainty_flag, filter_covariance_flag, smoother_redux, occbin_)
% Computes the diffuse kalman smoother in the case of a singular var-cov matrix.
% Univariate treatment of multivariate time series.
%
% INPUTS
%    a_initial:mm*1 vector of initial states
%    T:        mm*mm matrix     state transition matrix
%    Z:        pp*mm matrix     selector matrix for observables in augmented state vector
%    R:        mm*rr matrix     second matrix of the state equation relating the structural innovations to the state variables
%    Q:        rr*rr matrix     covariance matrix of structural errors
%    H:        pp*1             vector of variance of measurement errors
%    Pinf1:    mm*mm diagonal matrix with with q ones and m-q zeros
%    Pstar1:   mm*mm variance-covariance matrix with stationary variables
%    Y:        pp*1 vector
%    pp:       number of observed variables
%    mm:       number of state variables
%    smpl:     sample size
%    data_index:                [cell]      1*smpl cell of column vectors of indices.
%    nk:                        number of forecasting periods
%    kalman_tol:                tolerance for zero divider
%    diffuse_kalman_tol:        tolerance for zero divider
%    decomp_flag:               if true, compute filter decomposition
%    state_uncertainty_flag:    if true, compute uncertainty about smoothed
%                               state estimate
%    decomp_flag:               if true, compute filter decomposition
%    filter_covariance_flag:    if true, compute filter covariance
%    smoother_redux:            if true, compute smoother on restricted
%                               state space, recover static variables from this
%
% OUTPUTS
%    alphahat: smoothed state variables (a_{t|T})
%    epsilonhat: measurement errors
%    etahat:   smoothed shocks
%    a:        matrix of updated variables (a_{t|t})
%    aK:       3D array of k step ahead filtered state variables (a_{t+k|t})
%              (meaningless for periods 1:d)
%    P1:        3D array of one-step ahead forecast error variance
%              matrices
%    PK:       4D array of k-step ahead forecast error variance
%              matrices (meaningless for periods 1:d)
%    decomp:   decomposition of the effect of shocks on filtered values
%    V:        3D array of state uncertainty matrices
%    aalphahat:     filtered states in t-1|t
%    eetahat:       updated shocks in t|t
%    d:             number of diffuse periods
%
% Notes:
%   Outputs are stored in decision-rule order, i.e. to get variables in order of declaration
%   as in M_.endo_names, ones needs code along the lines of:
%   variables_declaration_order(dr.order_var,:) = alphahat
%
% Algorithm:
%
%   Uses the univariate filter as described in Durbin/Koopman (2012): "Time
%   Series Analysis by State Space Methods", Oxford University Press,
%   Second Edition, Ch. 6.4 + 7.2.5
%   and
%   Koopman/Durbin (2000): "Fast Filtering and Smoothing for Multivariatze State Space
%   Models", in Journal of Time Series Analysis, vol. 21(3), pp. 281-296.
%
% SPECIAL REQUIREMENTS
%   See "Filtering and Smoothing of State Vector for Diffuse State Space
%   Models", S.J. Koopman and J. Durbin (2003), in Journal of Time Series
%   Analysis, vol. 24(1), pp. 85-98.

% Copyright © 2004-2021 Dynare Team
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

% Modified by M. Ratto
% New output argument aK: 1-step to nk-stpe ahed predictions)
% New input argument nk: max order of predictions in aK

if size(H,2)>1
    error('missing_DiffuseKalmanSmootherH3_Z:: H is not a vector. This must not happens')
end

d = 0;
decomp = [];
spinf           = size(Pinf1);
spstar          = size(Pstar1);
v               = zeros(pp,smpl);
a               = zeros(mm,smpl);
a1              = zeros(mm,smpl+1);
a(:,1)          = a_initial;
a1(:,1)         = a_initial;
aK              = zeros(nk,mm,smpl+nk);

Fstar           = zeros(pp,smpl);
Finf            = zeros(pp,smpl);
Fi              = zeros(pp,smpl);
Ki              = zeros(mm,pp,smpl);
Kstar           = zeros(mm,pp,smpl);
Kinf            = zeros(spstar(1),pp,smpl);
P               = zeros(mm,mm,smpl+1);
P1              = P;
if filter_covariance_flag
    PK              = zeros(nk,mm,mm,smpl+nk);
else
    PK              = [];
end
Pstar           = zeros(spstar(1),spstar(2),smpl);
Pstar(:,:,1)    = Pstar1;
Pinf            = zeros(spinf(1),spinf(2),smpl);
Pinf(:,:,1)     = Pinf1;
Pstar1          = Pstar;
Pinf1           = Pinf;
rr              = size(Q,1); % number of structural shocks
isqvec = false;
if ndim(Q)>2
    Qvec = Q;
    Q=Q(:,:,1);
    isqvec = true;
end
QQ              = R*Q*transpose(R);
QRt             = Q*transpose(R);
alphahat        = zeros(mm,smpl);
etahat          = zeros(rr,smpl);
if smoother_redux
    aalphahat       = alphahat;
    eetahat         = etahat;
else
    aalphahat       = [];
    eetahat         = [];
end
epsilonhat      = zeros(rr,smpl);
r               = zeros(mm,smpl);
if state_uncertainty_flag
    if smoother_redux
        V               = zeros(mm+rr,mm+rr,smpl);
    else
        V               = zeros(mm,mm,smpl);
    end
    N               = zeros(mm,mm,smpl);
else
    V=[];
end

if ~occbin_.status
    isoccbin = 0;
    C=0;
    TT=[];
    RR=[];
    CC=[];
else
    isoccbin = 1;
    Qt = repmat(Q,[1 1 3]);
    options_=occbin_.info{1};
    oo_=occbin_.info{2};
    M_=occbin_.info{3};
    occbin_options=occbin_.info{4};
    opts_regime = occbin_options.opts_regime;
    %     first_period_occbin_update = inf;
    if isfield(opts_regime,'regime_history') && ~isempty(opts_regime.regime_history)
        opts_regime.regime_history=[opts_regime.regime_history(1) opts_regime.regime_history];
    else
        opts_regime.binding_indicator=zeros(smpl+2,M_.occbin.constraint_nbr);
    end
    occbin_options.opts_regime = opts_regime;
    [~, ~, ~, regimes_] = occbin.check_regimes([], [], [], opts_regime, M_, oo_, options_);
    if length(occbin_.info)>4
        if length(occbin_.info)==6 && options_.smoother_redux
            TT=repmat(T,1,1,smpl+1);
            RR=repmat(R,1,1,smpl+1);
            CC=repmat(zeros(mm,1),1,smpl+1);
            T0=occbin_.info{5};
            R0=occbin_.info{6};
        else
            
            TT=occbin_.info{5};
            RR=occbin_.info{6};
            CC=occbin_.info{7};
            %         TT = cat(3,TT,T);
            %         RR = cat(3,RR,R);
            %         CC = cat(2,CC,zeros(mm,1));
            if options_.smoother_redux
                my_order_var = oo_.dr.restrict_var_list;
                CC = CC(my_order_var,:);
                RR = RR(my_order_var,:,:);
                TT = TT(my_order_var,my_order_var,:);
                T0=occbin_.info{8};
                R0=occbin_.info{9};
            end
            if size(TT,3)<(smpl+1)
                TT=repmat(T,1,1,smpl+1);
                RR=repmat(R,1,1,smpl+1);
                CC=repmat(zeros(mm,1),1,smpl+1);
            end
        end
        
    else
        TT=repmat(T,1,1,smpl+1);
        RR=repmat(R,1,1,smpl+1);
        CC=repmat(zeros(mm,1),1,smpl+1);
    end
    if ~smoother_redux
        T0=T;
        R0=R;
        
    end
    if ~isinf(occbin_options.first_period_occbin_update)
        % initialize state matrices (otherwise they are set to 0 for
        % t<first_period_occbin_update!)
        TTT=repmat(T0,1,1,smpl+1);
        RRR=repmat(R0,1,1,smpl+1);
        CCC=repmat(zeros(length(T0),1),1,smpl+1);
    end
    
end

t = 0;
icc=0;
if ~isempty(Pinf(:,:,1))
    newRank = rank(Z*Pinf(:,:,1)*Z',diffuse_kalman_tol);
else
    newRank = rank(Pinf(:,:,1),diffuse_kalman_tol);
end
while newRank && t < smpl
    t = t+1;
    a(:,t) = a1(:,t);
    Pstar1(:,:,t) = Pstar(:,:,t);
    Pinf1(:,:,t) = Pinf(:,:,t);
    di = data_index{t}';
    for i=di
        Zi = Z(i,:);
        v(i,t)      = Y(i,t)-Zi*a(:,t);                                     % nu_{t,i} in 6.13 in DK (2012)
        Fstar(i,t)  = Zi*Pstar(:,:,t)*Zi' +H(i);                            % F_{*,t} in 5.7 in DK (2012), relies on H being diagonal
        Finf(i,t)   = Zi*Pinf(:,:,t)*Zi';                                   % F_{\infty,t} in 5.7 in DK (2012)
        Kstar(:,i,t) = Pstar(:,:,t)*Zi';                                    % KD (2000), eq. (15)
        if Finf(i,t) > diffuse_kalman_tol && newRank                        % F_{\infty,t,i} = 0, use upper part of bracket on p. 175 DK (2012) for w_{t,i}
            icc=icc+1;
            Kinf(:,i,t)       = Pinf(:,:,t)*Zi';                            % KD (2000), eq. (15)
            Kinf_Finf         = Kinf(:,i,t)/Finf(i,t);
            a(:,t)            = a(:,t) + Kinf_Finf*v(i,t);                  % KD (2000), eq. (16)
            Pstar(:,:,t)      = Pstar(:,:,t) + ...
                Kinf(:,i,t)*Kinf_Finf'*(Fstar(i,t)/Finf(i,t)) - ...
                Kstar(:,i,t)*Kinf_Finf' - ...
                Kinf_Finf*Kstar(:,i,t)';                                    % KD (2000), eq. (16)
            Pinf(:,:,t)       = Pinf(:,:,t) - Kinf(:,i,t)*Kinf(:,i,t)'/Finf(i,t); % KD (2000), eq. (16)
        elseif Fstar(i,t) > kalman_tol
            a(:,t)            = a(:,t) + Kstar(:,i,t)*v(i,t)/Fstar(i,t);    % KD (2000), eq. (17)
            Pstar(:,:,t)      = Pstar(:,:,t) - Kstar(:,i,t)*Kstar(:,i,t)'/Fstar(i,t);   % KD (2000), eq. (17)
            % Pinf is passed through unaltered, see eq. (17) of
            % Koopman/Durbin (2000)
        else
            % do nothing as a_{t,i+1}=a_{t,i} and P_{t,i+1}=P_{t,i}, see
            % p. 157, DK (2012)
        end
    end
    if newRank
        if ~isempty(Pinf(:,:,t))
            oldRank = rank(Z*Pinf(:,:,t)*Z',diffuse_kalman_tol);
        else
            oldRank = rank(Pinf(:,:,t),diffuse_kalman_tol);
        end
    else
        oldRank = 0;
    end
    if isoccbin,
        TT(:,:,t+1)=  T;
        RR(:,:,t+1)=  R;
    end
    a1(:,t+1) = T*a(:,t);
    aK(1,:,t+1) = a1(:,t+1);
    for jnk=2:nk
        aK(jnk,:,t+jnk) = T*dynare_squeeze(aK(jnk-1,:,t+jnk-1));
    end
    if isqvec
        QQ = R*Qvec(:,:,t+1)*transpose(R);
    end
    Pstar(:,:,t+1) = T*Pstar(:,:,t)*T'+ QQ;
    Pinf(:,:,t+1) = T*Pinf(:,:,t)*T';
    if newRank
        if ~isempty(Pinf(:,:,t+1))
            newRank = rank(Z*Pinf(:,:,t+1)*Z',diffuse_kalman_tol);
        else
            newRank = rank(Pinf(:,:,t+1),diffuse_kalman_tol);
        end
    end
    if oldRank ~= newRank
        disp('univariate_diffuse_kalman_filter:: T does influence the rank of Pinf!')
        disp('This may happen for models with order of integration >1.')
    end
end

if isoccbin
    first_period_occbin_update = max(t+2,occbin_options.first_period_occbin_update);
    if occbin_options.opts_regime.waitbar
        hh = dyn_waitbar(0,'Occbin: Piecewise Kalman Filter');
        set(hh,'Name','Occbin: Piecewise Kalman Filter.');
        waitbar_indicator=1;
    else
        waitbar_indicator=0;
    end
else
    first_period_occbin_update = inf;
    waitbar_indicator=0;
end
d = t;
P(:,:,d+1) = Pstar(:,:,d+1);
Fstar = Fstar(:,1:d);
Finf = Finf(:,1:d);
Kstar = Kstar(:,:,1:d);
Pstar = Pstar(:,:,1:d);
Pinf  = Pinf(:,:,1:d);
Pstar1 = Pstar1(:,:,1:d);
Pinf1  = Pinf1(:,:,1:d);
notsteady = 1;
while notsteady && t<smpl
    t = t+1;
    a(:,t) = a1(:,t);
    P1(:,:,t) = P(:,:,t);
    di = data_index{t}';
    if t>=first_period_occbin_update
        if waitbar_indicator
            dyn_waitbar(t/smpl, hh, sprintf('Period %u of %u', t,smpl));
        end
        if isqvec
            Qt = Qvec(:,:,t-1:t+1);
        end
        occbin_options.opts_regime.waitbar=0;
        [ax, a1x, Px, P1x, vx, Fix, Kix, Tx, Rx, Cx, tmp, error_flag, M_, aha, etaha,TTx,RRx,CCx] = occbin.kalman_update_algo_3(a(:,t-1),a1(:,t-1:t),P(:,:,t-1),P1(:,:,t-1:t),data_index(t-1:t),Z,v(:,t-1:t),Fi(:,t-1),Ki(:,:,t-1),Y(:,t-1:t),H,Qt,T0,R0,TT(:,:,t-1:t),RR(:,:,t-1:t),CC(:,t-1:t),regimes_(t:t+1),M_,oo_,options_,occbin_options,kalman_tol,nk);
        if ~error_flag
            regimes_(t:t+2)=tmp;
        else
            varargout{1} = [];
            varargout{2} = [];
            varargout{3} = [];
            varargout{4} = [];
            return
        end

        if smoother_redux
            aalphahat(:,t-1) = aha(:,1);
            eetahat(:,t) = etaha(:,2);
        end
        a(:,t) = ax(:,1);
        a1(:,t) = a1x(:,2);
        a1(:,t+1) = ax(:,2);
        v(di,t) = vx(di,2);
        Fi(di,t) = Fix(di,2);
        Ki(:,di,t) = Kix(:,di,2);
        TT(:,:,t:t+1) = Tx(:,:,1:2);
        RR(:,:,t:t+1) = Rx(:,:,1:2);
        CC(:,t:t+1) = Cx(:,1:2);
        TTT(:,:,t)=TTx;
        RRR(:,:,t)=RRx;
        CCC(:,t)=CCx;
        P(:,:,t) = Px(:,:,1);
        P1(:,:,t) = P1x(:,:,2);
        P(:,:,t+1) = Px(:,:,2);
        for jnk=1:nk
            PK(jnk,:,:,t+jnk) = Px(:,:,1+jnk);
            aK(jnk,:,t+jnk) = ax(:,1+jnk);
        end
    else
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
        if isqvec
            QQ = R*Qvec(:,:,t)*transpose(R);
        end
        if smoother_redux
            ri=zeros(mm,1);
            for st=t:-1:max(d+1,t-1)
                di = flipud(data_index{st})';
                for i = di
                    if Fi(i,st) > kalman_tol
                        Li = eye(mm)-Ki(:,i,st)*Z(i,:)/Fi(i,st);
                        ri = Z(i,:)'/Fi(i,st)*v(i,st)+Li'*ri;                             % DK (2012), 6.15, equation for r_{t,i-1}
                    end
                end
                if st==t-1
                    aalphahat(:,st) = a1(:,st) + P1(:,:,st)*ri;
                else
                    if isoccbin
                        if isqvec
                            QRt = Qvec(:,:,st)*transpose(RR(:,:,st));
                        else
                            QRt = Q*transpose(RR(:,:,st));
                        end
                        T = TT(:,:,st);
                    else
                        if isqvec
                            QRt = Qvec(:,:,st)*transpose(R);
                        end
                    end
                    eetahat(:,st) = QRt*ri;
                end
                ri = T'*ri;                                                             % KD (2003), eq. (23), equation for r_{t-1,p_{t-1}}
            end
        end
        if isoccbin
            if isqvec
                QQ = RR(:,:,t+1)*Qvec(:,:,t+1)*transpose(RR(:,:,t+1));
            else
                QQ = RR(:,:,t+1)*Q*transpose(RR(:,:,t+1));
            end
            T = TT(:,:,t+1);
            C = CC(:,t+1);
        else
            if isqvec
                QQ = R*Qvec(:,:,t+1)*transpose(R);
            end
        end
        a1(:,t+1) = T*a(:,t)+C;                                                 %transition according to (6.14) in DK (2012)
        P(:,:,t+1) = T*P(:,:,t)*T' + QQ;                                        %transition according to (6.14) in DK (2012)
        if filter_covariance_flag
            Pf          = P(:,:,t+1);
        end
        aK(1,:,t+1) = a1(:,t+1);
        if ~isempty(nk) && nk>1 && isoccbin && (t>=first_period_occbin_update || isinf(first_period_occbin_update))
            opts_simul = occbin_options.opts_regime;
            opts_simul.SHOCKS = zeros(nk,M_.exo_nbr);
            if smoother_redux
                tmp=zeros(M_.endo_nbr,1);
                tmp(oo_.dr.restrict_var_list)=a(:,t);
                opts_simul.endo_init = tmp(oo_.dr.inv_order_var);
            else
                opts_simul.endo_init = a(oo_.dr.inv_order_var,t);
            end
            opts_simul.init_regime = []; %regimes_(t);
            opts_simul.waitbar=0;
            options_.occbin.simul=opts_simul;
            [~, out, ss] = occbin.solver(M_,oo_,options_);
        end
        for jnk=1:nk
            if filter_covariance_flag
                if jnk>1
                    Pf = T*Pf*T' + QQ;
                end
                PK(jnk,:,:,t+jnk) = Pf;
            end
            if isoccbin && (t>=first_period_occbin_update || isinf(first_period_occbin_update))
                if smoother_redux
                    aK(jnk,:,t+jnk) = out.piecewise(jnk,oo_.dr.order_var(oo_.dr.restrict_var_list)) - out.ys(oo_.dr.order_var(oo_.dr.restrict_var_list))';
                else
                    aK(jnk,oo_.dr.inv_order_var,t+jnk) = out.piecewise(jnk,:) - out.ys';
                end
            elseif jnk>1
                aK(jnk,:,t+jnk) = T*dynare_squeeze(aK(jnk-1,:,t+jnk-1));
            end
        end
    end
end
if waitbar_indicator
    dyn_waitbar_close(hh); 
end

P1(:,:,t+1) = P(:,:,t+1);

if ~isinf(first_period_occbin_update) && isoccbin
    regimes_ = regimes_(2:smpl+1);
else
    regimes_ = struct();
    TTT=TT;
    RRR=RR;
    CCC=CC;
    %     return
end
varargout{1} = regimes_;
varargout{2} = TTT;
varargout{3} = RRR;
varargout{4} = CCC;
% $$$ P_s=tril(P(:,:,t))+tril(P(:,:,t),-1)';
% $$$ P1_s=tril(P1(:,:,t))+tril(P1(:,:,t),-1)';
% $$$ Fi_s = Fi(:,t);
% $$$ Ki_s = Ki(:,:,t);
% $$$ L_s  =Li(:,:,:,t);
% $$$ if t<smpl
% $$$   P  = cat(3,P(:,:,1:t),repmat(P_s,[1 1 smpl-t]));
% $$$   P1  = cat(3,P1(:,:,1:t),repmat(P1_s,[1 1 smpl-t]));
% $$$   Fi = cat(2,Fi(:,1:t),repmat(Fi_s,[1 1 smpl-t]));
% $$$   Li  = cat(4,Li(:,:,:,1:t),repmat(L_s,[1 1 smpl-t]));
% $$$   Ki  = cat(3,Ki(:,:,1:t),repmat(Ki_s,[1 1 smpl-t]));
% $$$ end
% $$$ while t<smpl
% $$$   t=t+1;
% $$$   a(:,t) = a1(:,t);
% $$$   di = data_index{t}';
% $$$   for i=di
% $$$     Zi = Z(i,:);
% $$$     v(i,t)      = Y(i,t) - Zi*a(:,t);
% $$$     if Fi_s(i) > kalman_tol
% $$$       a(:,t) = a(:,t) + Ki_s(:,i)*v(i,t)/Fi_s(i);
% $$$     end
% $$$   end
% $$$   a1(:,t+1) = T*a(:,t);
% $$$   Pf          = P(:,:,t);
% $$$   for jnk=1:nk,
% $$$       Pf = T*Pf*T' + QQ;
% $$$       aK(jnk,:,t+jnk) = T^jnk*a(:,t);
% $$$       PK(jnk,:,:,t+jnk) = Pf;
% $$$   end
% $$$ end

%% do backward pass
ri=zeros(mm,1);
if state_uncertainty_flag
    Ni=zeros(mm,mm);
end
t = smpl+1;
while t > d+1
    t = t-1;
    di = flipud(data_index{t})';
    for i = di
        if Fi(i,t) > kalman_tol
            Li = eye(mm)-Ki(:,i,t)*Z(i,:)/Fi(i,t);
            ri = Z(i,:)'/Fi(i,t)*v(i,t)+Li'*ri;                             % DK (2012), 6.15, equation for r_{t,i-1}
            if state_uncertainty_flag
                Ni = Z(i,:)'/Fi(i,t)*Z(i,:)+Li'*Ni*Li;                      % KD (2000), eq. (23)
            end
        end
    end
    r(:,t) = ri;                                                            % DK (2012), below 6.15, r_{t-1}=r_{t,0}
    alphahat(:,t) = a1(:,t) + P1(:,:,t)*r(:,t);
    if isoccbin
        if isqvec
            QRt = Qvec(:,:,t)*transpose(RR(:,:,t));
        else
            QRt = Q*transpose(RR(:,:,t));
        end
        R = RR(:,:,t);
        T = TT(:,:,t);
    else
        if isqvec
            QRt             = Qvec(:,:,t)*transpose(R);
        end
    end
    etahat(:,t) = QRt*r(:,t);
    ri = T'*ri;                                                             % KD (2003), eq. (23), equation for r_{t-1,p_{t-1}}
    if state_uncertainty_flag
        N(:,:,t) = Ni;                                                          % DK (2012), below 6.15, N_{t-1}=N_{t,0}
        if smoother_redux
            ptmp = [P1(:,:,t) R*Q; (R*Q)' Q];
            ntmp = [N(:,:,t) zeros(mm,rr); zeros(rr,mm+rr)];
            V(:,:,t)    = ptmp - ptmp*ntmp*ptmp;
        else
            V(:,:,t) = P1(:,:,t)-P1(:,:,t)*N(:,:,t)*P1(:,:,t);                      % KD (2000), eq. (7) with N_{t-1} stored in N(:,:,t)
        end
        Ni = T'*Ni*T;                                                           % KD (2000), eq. (23), equation for N_{t-1,p_{t-1}}
    end
end
if d
    r0 = zeros(mm,d);
    r0(:,d) = ri;
    r1 = zeros(mm,d);
    if state_uncertainty_flag
        %N_0 at (d+1) is N(d+1), so we can use N for continuing and storing N_0-recursion
        N_0=zeros(mm,mm,d);   %set N_1_{d}=0, below  KD (2000), eq. (24)
        N_0(:,:,d) = Ni;
        N_1=zeros(mm,mm,d);   %set N_1_{d}=0, below  KD (2000), eq. (24)
        N_2=zeros(mm,mm,d);   %set N_2_{d}=0, below  KD (2000), eq. (24)
    end
    for t = d:-1:1
        di = flipud(data_index{t})';
        for i = di
            if Finf(i,t) > diffuse_kalman_tol
                % recursions need to be from highest to lowest term in order to not
                % overwrite lower terms still needed in this step
                Linf    = eye(mm) - Kinf(:,i,t)*Z(i,:)/Finf(i,t);
                L0      = (Kinf(:,i,t)*(Fstar(i,t)/Finf(i,t))-Kstar(:,i,t))*Z(i,:)/Finf(i,t);
                r1(:,t) = Z(i,:)'*v(i,t)/Finf(i,t) + ...
                    L0'*r0(:,t) + ...
                    Linf'*r1(:,t);   % KD (2000), eq. (25) for r_1
                r0(:,t) = Linf'*r0(:,t);   % KD (2000), eq. (25) for r_0
                if state_uncertainty_flag
                    N_2(:,:,t)=Z(i,:)'/Finf(i,t)^2*Z(i,:)*Fstar(i,t) ...
                        + Linf'*N_2(:,:,t)*Linf...
                        + Linf'*N_1(:,:,t)*L0...
                        + L0'*N_1(:,:,t)'*Linf...
                        + L0'*N_0(:,:,t)*L0;                                    % DK (2012), eq. 5.29
                    N_1(:,:,t)=Z(i,:)'/Finf(i,t)*Z(i,:)+Linf'*N_1(:,:,t)*Linf...
                        +L0'*N_0(:,:,t)*Linf;                                   % DK (2012), eq. 5.29; note that, compared to DK (2003) this drops the term (L_1'*N(:,:,t+1)*Linf(:,:,t))' in the recursion due to it entering premultiplied by Pinf when computing V, and Pinf*Linf'*N=0
                    N_0(:,:,t)=Linf'*N_0(:,:,t)*Linf;                           % DK (2012), eq. 5.19, noting that L^(0) is named Linf
                end
            elseif Fstar(i,t) > kalman_tol % step needed whe Finf == 0
                L_i=eye(mm) - Kstar(:,i,t)*Z(i,:)/Fstar(i,t);
                r0(:,t) = Z(i,:)'/Fstar(i,t)*v(i,t)+L_i'*r0(:,t);           % propagate r0 and keep r1 fixed
                if state_uncertainty_flag
                    N_0(:,:,t)=Z(i,:)'/Fstar(i,t)*Z(i,:)+L_i'*N_0(:,:,t)*L_i;   % propagate N_0 and keep N_1 and N_2 fixed
                end
            end
        end
        alphahat(:,t) = a1(:,t) + Pstar1(:,:,t)*r0(:,t) + Pinf1(:,:,t)*r1(:,t); % KD (2000), eq. (26)
        r(:,t)        = r0(:,t);
        if isoccbin
            if isqvec
                QRt = Qvec(:,:,t)*transpose(RR(:,:,t));
            else
                QRt = Q*transpose(RR(:,:,t));
            end
            R = RR(:,:,t);
            T = TT(:,:,t);
        else
            if isqvec
                QRt             = Qvec(:,:,t)*transpose(R);
            end
        end
        etahat(:,t)   = QRt*r(:,t);                                         % KD (2000), eq. (27)
        if state_uncertainty_flag
            if smoother_redux
                pstmp = [Pstar(:,:,t) R*Q; (R*Q)' Q];
                pitmp = [Pinf(:,:,t) zeros(mm,rr); zeros(rr,mm+rr)];
                ntmp0 = [N_0(:,:,t) zeros(mm,rr); zeros(rr,mm+rr)];
                ntmp1 = [N_1(:,:,t) zeros(mm,rr); zeros(rr,mm+rr)];
                ntmp2 = [N_2(:,:,t) zeros(mm,rr); zeros(rr,mm+rr)];
                V(:,:,t)    = pstmp - pstmp*ntmp0*pstmp...
                    -(pitmp*ntmp1*pstmp)'...
                    - pitmp*ntmp1*pstmp...
                    - pitmp*ntmp2*Pinf(:,:,t);                                   % DK (2012), eq. 5.30
                
            else
                V(:,:,t)=Pstar(:,:,t)-Pstar(:,:,t)*N_0(:,:,t)*Pstar(:,:,t)...
                    -(Pinf(:,:,t)*N_1(:,:,t)*Pstar(:,:,t))'...
                    - Pinf(:,:,t)*N_1(:,:,t)*Pstar(:,:,t)...
                    - Pinf(:,:,t)*N_2(:,:,t)*Pinf(:,:,t);                       % DK (2012), eq. 5.30
            end
        end
        if t > 1
            r0(:,t-1) = T'*r0(:,t);                                         % KD (2000), below eq. (25) r_{t-1,p_{t-1}}=T'*r_{t,0}
            r1(:,t-1) = T'*r1(:,t);                                         % KD (2000), below eq. (25) r_{t-1,p_{t-1}}=T'*r_{t,0}
            if state_uncertainty_flag
                N_0(:,:,t-1)= T'*N_0(:,t)*T;                                % KD (2000), below eq. (25) N_{t-1,p_{t-1}}=T'*N_{t,0}*T
                N_1(:,:,t-1)= T'*N_1(:,t)*T;                                % KD (2000), below eq. (25) N^1_{t-1,p_{t-1}}=T'*N^1_{t,0}*T
                N_2(:,:,t-1)= T'*N_2(:,t)*T;                                % KD (2000), below eq. (25) N^2_{t-1,p_{t-1}}=T'*N^2_{t,0}*T
            end
        end
    end
else
    alphahat0 = 0*a1(:,1) + P1(:,:,1)*ri;
end

if decomp_flag
    decomp = zeros(nk,mm,rr,smpl+nk);
    ZRQinv = inv(Z*QQ*Z');
    for t = max(d,1):smpl
        ri_d = zeros(mm,1);
        di = flipud(data_index{t})';
        for i = di
            if Fi(i,t) > kalman_tol
                ri_d = Z(i,:)'/Fi(i,t)*v(i,t)+ri_d-Ki(:,i,t)'*ri_d/Fi(i,t)*Z(i,:)';
            end
        end
        
        % calculate eta_tm1t
        if isoccbin
            if isqvec
                QRt = Qvec(:,:,t)*transpose(RR(:,:,t));
            else
                QRt = Q*transpose(RR(:,:,t));
            end
            R = RR(:,:,t);
            T = TT(:,:,t);
        else
            if isqvec
                QRt = Qvec(:,:,t)*transpose(R);
            end
        end
        eta_tm1t = QRt*ri_d;
        % calculate decomposition
        Ttok = eye(mm,mm);
        AAA = P1(:,:,t)*Z'*ZRQinv*Z*R;
        for h = 1:nk
            BBB = Ttok*AAA;
            for j=1:rr
                decomp(h,:,j,t+h) = eta_tm1t(j)*BBB(:,j);
            end
            Ttok = T*Ttok;
        end
    end
end

epsilonhat = Y - Z*alphahat;


if (d==smpl)
    warning(['missing_DiffuseKalmanSmootherH3_Z:: There isn''t enough information to estimate the initial conditions of the nonstationary variables']);
    return
end
