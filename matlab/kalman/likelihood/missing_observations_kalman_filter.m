function  [LIK, lik, a, P] = missing_observations_kalman_filter(data_index,number_of_observations,no_more_missing_observations,Y,start,last,a,P,kalman_tol,riccati_tol,rescale_prediction_error_covariance,presample,T,Q,R,H,Z,mm,pp,rr,Zflag,diffuse_periods,occbin_)
% [LIK, lik, a, P] = missing_observations_kalman_filter(data_index,number_of_observations,no_more_missing_observations,Y,start,last,a,P,kalman_tol,riccati_tol,rescale_prediction_error_covariance,presample,T,Q,R,H,Z,mm,pp,rr,Zflag,diffuse_periods,occbin_)
% Computes the likelihood of a state space model in the case with missing observations.
%
% INPUTS
%    data_index                   [cell]      1*smpl cell of column vectors of indices.
%    number_of_observations       [integer]   scalar.
%    no_more_missing_observations [integer]   scalar.
%    Y                            [double]    pp*smpl matrix of data.
%    start                        [integer]   scalar, index of the first observation.
%    last                         [integer]   scalar, index of the last observation.
%    a                            [double]    pp*1 vector, levels of the predicted initial state variables (E_{0}(alpha_1)).
%    P                            [double]    pp*pp matrix, covariance matrix of the initial state vector.
%    kalman_tol                   [double]    scalar, tolerance parameter (rcond).
%    riccati_tol                  [double]    scalar, tolerance parameter (riccati iteration).
%    presample                    [integer]   scalar, presampling if strictly positive.
%    T                            [double]    mm*mm transition matrix of the state equation.
%    Q                            [double]    rr*rr covariance matrix of the structural innovations.
%    R                            [double]    mm*rr matrix, mapping structural innovations to state variables.
%    H                            [double]    pp*pp (or 1*1 =0 if no measurement error) covariance matrix of the measurement errors.
%    Z                            [integer]   pp*1 vector of indices for the observed variables.
%    mm                           [integer]   scalar, dimension of the state vector.
%    pp                           [integer]   scalar, number of observed variables.
%    rr                           [integer]   scalar, number of structural innovations.
%
% OUTPUTS
%    LIK        [double]    scalar, MINUS loglikelihood
%    lik        [double]    vector, density of observations in each period.
%    a          [double]    mm*1 vector, current estimate of the state vector tomorrow (E_{T}(alpha_{T+1})).
%    P          [double]    mm*mm matrix, covariance matrix of the states.
%
%
% NOTES
%   The vector "lik" is used to evaluate the jacobian of the likelihood.

% Copyright © 2004-2023 Dynare Team
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

% Set defaults
if nargin<20
    Zflag = 0;
    diffuse_periods = 0;
end

if nargin<21
    diffuse_periods = 0;
end

if isempty(Zflag)
    Zflag = 0;
end

if isempty(diffuse_periods)
    diffuse_periods = 0;
end

if isequal(H,0)
    H = zeros(pp,pp);
end

% Get sample size.
smpl = last-start+1;

% Initialize some variables.
dF   = 1;
isqvec = false;
if ndims(Q)>2
    Qvec = Q;
    Q=Q(:,:,1);
    isqvec = true;
end
QQ   = R*Q*transpose(R);   % Variance of R times the vector of structural innovations.
t    = start;              % Initialization of the time index.
lik  = zeros(smpl,1);      % Initialization of the vector gathering the densities.
LIK  = Inf;                % Default value of the log likelihood.
oldK = Inf;
notsteady   = 1;
F_singular  = true;
s = 0;
rescale_prediction_error_covariance0=rescale_prediction_error_covariance;
if occbin_.status
    base_regime = struct();
    Qt = repmat(Q,[1 1 3]);
    a0 = zeros(mm,last);
    a1 = zeros(mm,last);
    P0 = zeros(mm,mm,last);
    P1 = zeros(mm,mm,last);
    vv = zeros(pp,last);

    options_=occbin_.info{1};
    dr=occbin_.info{2};
    endo_steady_state=occbin_.info{3};
    exo_steady_state=occbin_.info{4};
    exo_det_steady_state=occbin_.info{5};
    M_=occbin_.info{6};
    occbin_options=occbin_.info{7};
    occbin_options.opts_simul.SHOCKS = [];
    opts_regime.regime_history = occbin_options.opts_regime.init_regime_history;
    opts_regime.binding_indicator = occbin_options.opts_regime.init_binding_indicator;
    if t>1
        first_period_occbin_update = max(t+1,options_.occbin.likelihood.first_period_occbin_update);
    else
        first_period_occbin_update = options_.occbin.likelihood.first_period_occbin_update;
    end
    if isempty(opts_regime.binding_indicator) && isempty(opts_regime.regime_history)
        opts_regime.binding_indicator=zeros(last+2,M_.occbin.constraint_nbr);
    end
    [~, ~, ~, regimes_] = occbin.check_regimes([], [], [], opts_regime, M_, options_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
    if length(occbin_.info)>7
        TT=occbin_.info{8};
        RR=occbin_.info{9};
        CC=occbin_.info{10};
        T0=occbin_.info{11};
        R0=occbin_.info{12};
        TT = cat(3,TT,T);
        RR = cat(3,RR,R);
        CC = cat(2,CC,zeros(mm,1));
        if size(TT,3)<(last+1)
            TT=repmat(T,1,1,last+1);
            RR=repmat(R,1,1,last+1);
            CC=repmat(zeros(mm,1),1,last+1);
        end

    end
    if M_.occbin.constraint_nbr==1
        base_regime.regime = 0;
        base_regime.regimestart = 1;
    else
        base_regime.regime1 = 0;
        base_regime.regimestart1 = 1;
        base_regime.regime2 = 0;
        base_regime.regimestart2 = 1;
    end
else
    first_period_occbin_update = inf;
    C=0;
end

while notsteady && t<=last
    if occbin_.status
        a1(:,t) = a;
        P1(:,:,t) = P;
        C = CC(:,t+1);
        R = RR(:,:,t+1);
        T = TT(:,:,t+1);
        if ~(isqvec)
            QQ = R*Q*transpose(R);   % Variance of R times the vector of structural innovations.
        end
        if t==1
            Pinit = P1(:,:,1);
            ainit = a1(:,1);
        end
    end
    s  = t-start+1;
    d_index = data_index{t};
    if isqvec
        QQ = R*Qvec(:,:,t+1)*transpose(R);
    end
    if isempty(d_index)
        a = T*a;
        P = T*P*transpose(T)+QQ;
    else
        % Compute the prediction error and its variance
        if Zflag
            z = Z(d_index,:);
            v = Y(d_index,t)-z*a;
            F = z*P*z' + H(d_index,d_index);
        else
            z = Z(d_index);
            v = Y(d_index,t) - a(z);
            F = P(z,z) + H(d_index,d_index);
        end
        badly_conditioned_F = false;
        if rescale_prediction_error_covariance
            sig=sqrt(diag(F));
            if any(diag(F)<kalman_tol) || rcond(F./(sig*sig'))<kalman_tol
                badly_conditioned_F = true;
            end
        else
            if rcond(F)<kalman_tol
                sig=sqrt(diag(F));
                if any(diag(F)<kalman_tol) || rcond(F./(sig*sig'))<kalman_tol
                    badly_conditioned_F = true;
                else
                    rescale_prediction_error_covariance=1;
                end
                %                 badly_conditioned_F = true;
            end
        end
        if badly_conditioned_F && (~occbin_.status || (occbin_.status && t<first_period_occbin_update))
            % if ~all(abs(F(:))<kalman_tol), then use univariate filter, otherwise this is a
            % pathological case and the draw is discarded
            return
        else
            F_singular = false;
        end
        if  ~occbin_.status || (occbin_.status && (options_.occbin.likelihood.use_updated_regime==0 || t<first_period_occbin_update))
            if rescale_prediction_error_covariance
                log_dF = log(det(F./(sig*sig')))+2*sum(log(sig));
                iF = inv(F./(sig*sig'))./(sig*sig');
                rescale_prediction_error_covariance=rescale_prediction_error_covariance0;
            else
                log_dF = log(det(F));
                iF = inv(F);
            end
            lik(s) = log_dF + transpose(v)*iF*v + length(d_index)*log(2*pi);
            if t<first_period_occbin_update
                if Zflag
                    K = P*z'*iF;
                    if occbin_.status
                        P0(:,:,t) = (P-K*z*P);
                    end

                    P = T*(P-K*z*P)*transpose(T)+QQ;
                else
                    K = P(:,z)*iF;
                    if occbin_.status
                        P0(:,:,t) = (P-K*P(z,:));
                    end
                    P = T*(P-K*P(z,:))*transpose(T)+QQ;
                end
                if occbin_.status
                    a0(:,t) = (a+K*v);
                    vv(d_index,t) = v;
                end
                a = T*(a+K*v)+C;
                if t>=no_more_missing_observations && ~isqvec && ~occbin_.status
                    notsteady = max(abs(K(:)-oldK))>riccati_tol;
                    oldK = K(:);
                end
            end
        end
    end
    if occbin_.status && t>=first_period_occbin_update

        occbin_options.opts_simul.waitbar=0;
        if t==1
            if isqvec
                Qt = cat(3,Q,Qvec(:,:,t:t+1));
            end
            a00 = ainit;
            a10 = [a00 a(:,1)];
            P00 = Pinit;
            P10 = P1(:,:,[1 1]);
            data_index0{1}=[];
            data_index0(2)=data_index(1);
            v0(:,2)=vv(:,1);
            Y0(:,2)=Y(:,1);
            Y0(:,1)=nan;
            TT01 = cat(3,T,TT(:,:,1));
            RR01 = cat(3,R,RR(:,:,1));
            CC01 = zeros(size(CC,1),2);
            CC01(:,2) = CC(:,1);
            % insert here kalman update engine
            [ax, a1x, Px, P1x, vx, Tx, Rx, Cx, regx, info, M_, likx] = occbin.kalman_update_engine(a00, a10, P00, P10, t, data_index0, Z, v0, Y0, H, Qt, T0, R0, TT01, RR01, CC01, regimes_(t:t+1), base_regime, d_index, M_, dr, endo_steady_state,exo_steady_state,exo_det_steady_state, options_, occbin_options);
%            [ax, a1x, Px, P1x, vx, Tx, Rx, Cx, regimes_(t:t+2), info, M_, likx] = occbin.kalman_update_algo_1(a00, a10, P00, P10, data_index0, Z, v0, Y0, H, Qt, T0, R0, TT01, RR01, CC01, regimes_(t:t+1), M_, dr, endo_steady_state,exo_steady_state,exo_det_steady_state, options_, occbin_options);
        else
            if isqvec
                Qt = Qvec(:,:,t-1:t+1);
            end
            % insert here kalman update engine
            [ax, a1x, Px, P1x, vx, Tx, Rx, Cx, regx, info, M_, likx] = occbin.kalman_update_engine(a0(:,t-1),a1(:,t-1:t),P0(:,:,t-1),P1(:,:,t-1:t),t,data_index(t-1:t),Z,vv(:,t-1:t),Y(:,t-1:t),H,Qt,T0,R0,TT(:,:,t-1:t),RR(:,:,t-1:t),CC(:,t-1:t),regimes_(t:t+1),base_regime,d_index,M_,dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
%            [ax, a1x, Px, P1x, vx, Tx, Rx, Cx, regimes_(t:t+2), info, M_, likx] = occbin.kalman_update_algo_1(a0(:,t-1),a1(:,t-1:t),P0(:,:,t-1),P1(:,:,t-1:t),data_index(t-1:t),Z,vv(:,t-1:t),Y(:,t-1:t),H,Qt,T0,R0,TT(:,:,t-1:t),RR(:,:,t-1:t),CC(:,t-1:t),regimes_(t:t+1),M_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
        end
        if info
            if options_.debug
                fprintf('\nmissing_observations_kalman_filter:PKF failed in period %u with: %s\n', t, get_error_message(info,options_));
            end
            return
        end
        if options_.occbin.likelihood.use_updated_regime
            lik(s) = likx;
        end
        regimes_(t:t+2) = regx;
        a0(:,t) = ax(:,1);
        a1(:,t) = a1x(:,2);
        a = ax(:,2);
        vv(d_index,t) = vx(d_index,2);
        TT(:,:,t:t+1) = Tx;
        RR(:,:,t:t+1) = Rx;
        CC(:,t:t+1) = Cx;
        P0(:,:,t) = Px(:,:,1);
        P1(:,:,t) = P1x(:,:,2);
        P = Px(:,:,2);

    end
    t = t+1;
end

if F_singular
    error('The variance of the forecast error remains singular until the end of the sample')
end

% Divide by two.
lik(1:s) = .5*lik(1:s);

% Call steady state Kalman filter if needed.
if t<=last
    [~, lik(s+1:end)] = kalman_filter_ss(Y, t, last, a, T, K, iF, log_dF, Z, pp, Zflag);
end

% Compute minus the log-likelihood.
if presample>=diffuse_periods
    LIK = sum(lik(1+presample-diffuse_periods:end));
else
    LIK = sum(lik);
end
