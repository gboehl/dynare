function [AHess, DLIK, LIK] = AHessian(T,R,Q,H,P,Y,DT,DYss,DOm,DH,DP,start,mf,kalman_tol,riccati_tol)
% function [AHess, DLIK, LIK] = AHessian(T,R,Q,H,P,Y,DT,DYss,DOm,DH,DP,start,mf,kalman_tol,riccati_tol)
%
% computes the asymptotic hessian matrix of the log-likelihood function of
% a state space model (notation as in kalman_filter.m in DYNARE
% Thanks to  Nikolai Iskrev
%
% NOTE: the derivative matrices (DT,DR ...) are 3-dim. arrays with last
% dimension equal to the number of structural parameters

% Copyright © 2011-2017 Dynare Team
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


k = size(DT,3);                                 % number of structural parameters
smpl = size(Y,2);                               % Sample size.
pp   = size(Y,1);                               % Maximum number of observed variables.
mm   = size(T,2);                               % Number of state variables.
a    = zeros(mm,1);                             % State vector.
Om   = R*Q*transpose(R);                        % Variance of R times the vector of structural innovations.
t    = 0;                                       % Initialization of the time index.
oldK = 0;
notsteady   = 1;                                % Steady state flag.
F_singular  = 1;

lik  = zeros(smpl,1);                           % Initialization of the vector gathering the densities.
LIK  = Inf;                                     % Default value of the log likelihood.
if nargout > 1
    DLIK  = zeros(k,1);                             % Initialization of the score.
end
AHess  = zeros(k,k);                             % Initialization of the Hessian
Da    = zeros(mm,k);                             % State vector.
Dv = zeros(length(mf),k);

%     for ii = 1:k
%         DOm = DR(:,:,ii)*Q*transpose(R) + R*DQ(:,:,ii)*transpose(R) + R*Q*transpose(DR(:,:,ii));
%     end

while notsteady && t<smpl
    t  = t+1;
    v  = Y(:,t)-a(mf);
    F  = P(mf,mf) + H;
    if rcond(F) < kalman_tol
        if ~all(abs(F(:))<kalman_tol)
            return
        else
            a = T*a;
            P = T*P*transpose(T)+Om;
        end
    else
        F_singular = 0;
        iF     = inv(F);
        K      = P(:,mf)*iF;
        lik(t) = log(det(F))+transpose(v)*iF*v;

        [DK,DF,DP1] = computeDKalman(T,DT,DOm,P,DP,DH,mf,iF,K);

        for ii = 1:k
            Dv(:,ii)   = -Da(mf,ii) - DYss(mf,ii);
            Da(:,ii)   = DT(:,:,ii)*(a+K*v) + T*(Da(:,ii)+DK(:,:,ii)*v + K*Dv(:,ii));
            if t>=start && nargout > 1
                DLIK(ii,1)  = DLIK(ii,1) + trace( iF*DF(:,:,ii) ) + 2*Dv(:,ii)'*iF*v - v'*(iF*DF(:,:,ii)*iF)*v;
            end
        end
        vecDPmf = reshape(DP(mf,mf,:),[],k);
        %                     iPmf = inv(P(mf,mf));
        if t>=start
            AHess = AHess + Dv'*iF*Dv + .5*(vecDPmf' * kron(iF,iF) * vecDPmf);
        end
        a      = T*(a+K*v);
        P      = T*(P-K*P(mf,:))*transpose(T)+Om;
        DP     = DP1;
    end
    notsteady = max(max(abs(K-oldK))) > riccati_tol;
    oldK = K;
end

if F_singular
    error('The variance of the forecast error remains singular until the end of the sample')
end


if t < smpl
    t0 = t+1;
    while t < smpl
        t = t+1;
        v = Y(:,t)-a(mf);
        for ii = 1:k
            Dv(:,ii)   = -Da(mf,ii)-DYss(mf,ii);
            Da(:,ii)   = DT(:,:,ii)*(a+K*v) + T*(Da(:,ii)+DK(:,:,ii)*v + K*Dv(:,ii));
            if t>=start && nargout >1
                DLIK(ii,1)  = DLIK(ii,1) + trace( iF*DF(:,:,ii) ) + 2*Dv(:,ii)'*iF*v - v'*(iF*DF(:,:,ii)*iF)*v;
            end
        end
        if t>=start
            AHess = AHess + Dv'*iF*Dv;
        end
        a = T*(a+K*v);
        lik(t) = transpose(v)*iF*v;
    end
    AHess = AHess + .5*(smpl-t0+1)*(vecDPmf' * kron(iF,iF) * vecDPmf);
    if nargout > 1
        for ii = 1:k
            %             DLIK(ii,1)  = DLIK(ii,1) + (smpl-t0+1)*trace( iF*DF(:,:,ii) );
        end
    end
    lik(t0:smpl) = lik(t0:smpl) + log(det(F));
    %         for ii = 1:k;
    %             for jj = 1:ii
    %              H(ii,jj) = trace(iPmf*(.5*DP(mf,mf,ii)*iPmf*DP(mf,mf,jj) + Dv(:,ii)*Dv(:,jj)'));
    %             end
    %         end
end

AHess = -AHess;
if nargout > 1
    DLIK = DLIK/2;
end
% adding log-likelihhod constants
lik = (lik + pp*log(2*pi))/2;

LIK = sum(lik(start:end)); % Minus the log-likelihood.
                           % end of main function

function [DK,DF,DP1] = computeDKalman(T,DT,DOm,P,DP,DH,mf,iF,K)

k      = size(DT,3);
tmp    = P-K*P(mf,:);

for ii = 1:k
    DF(:,:,ii)  = DP(mf,mf,ii) + DH(:,:,ii);
    DiF(:,:,ii) = -iF*DF(:,:,ii)*iF;
    DK(:,:,ii)  = DP(:,mf,ii)*iF + P(:,mf)*DiF(:,:,ii);
    Dtmp        = DP(:,:,ii) - DK(:,:,ii)*P(mf,:) - K*DP(mf,:,ii);
    DP1(:,:,ii) = DT(:,:,ii)*tmp*T' + T*Dtmp*T' + T*tmp*DT(:,:,ii)' + DOm(:,:,ii);
end

% end of computeDKalman
