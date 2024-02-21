function [dLIK, dlikk, a, Pstar, llik] = univariate_kalman_filter_d(data_index, number_of_observations, no_more_missing_observations, Y, start, last, a, Pinf, Pstar, kalman_tol, diffuse_kalman_tol, riccati_tol, presample, T, R, Q, H, Z, mm, pp, rr)
% Computes the likelihood of a state space model (initialization with diffuse steps, univariate approach).

%@info:
%! @deftypefn {Function File} {[@var{LIK},@var{likk},@var{a},@var{Pstar}, @var{llik} ] =} univariate_kalman_filter_d (@var{data_index}, @var{number_of_observations},@var{no_more_missing_observations}, @var{Y}, @var{start}, @var{last}, @var{a}, @var{Pinf}, @var{Pstar}, @var{kalman_tol}, @var{riccati_tol},@var{presample},@var{T},@var{R},@var{Q},@var{H},@var{Z},@var{mm},@var{pp},@var{rr})
%! @anchor{univariate_kalman_filter_d}
%! @sp 1
%! Computes the likelihood of a state space model (initialization with diffuse steps, univariate approach).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item data_index
%! Matlab's cell, 1*T cell of column vectors of indices (in the vector of observed variables).
%! @item number_of_observations
%! Integer scalar, effective number of observations.
%! @item no_more_missing_observations
%! Integer scalar, date after which there is no more missing observation (it is then possible to switch to the steady state kalman filter).
%! @item Y
%! Matrix (@var{pp}*T) of doubles, data.
%! @item start
%! Integer scalar, first period.
%! @item last
%! Integer scalar, last period (@var{last}-@var{first} has to be inferior to T).
%! @item a
%! Vector (@var{mm}*1) of doubles, initial mean of the state vector.
%! @item Pinf
%! Matrix (@var{mm}*@var{mm}) of doubles, initial covariance matrix of the state vector (non stationary part).
%! @item Pstar
%! Matrix (@var{mm}*@var{mm}) of doubles, initial covariance matrix of the state vector (stationary part).
%! @item kalman_tol
%! Double scalar, tolerance parameter (rcond, inversibility of the covariance matrix of the prediction errors).
%! @item riccati_tol
%! Double scalar, tolerance parameter (iteration over the Riccati equation).
%! @item presample
%! Integer scalar, presampling if strictly positive (number of initial iterations to be discarded when evaluating the likelihood).
%! @item T
%! Matrix (@var{mm}*@var{mm}) of doubles, transition matrix of the state equation.
%! @item R
%! Matrix (@var{mm}*@var{rr}) of doubles, second matrix of the state equation relating the structural innovations to the state variables.
%! @item Q
%! Matrix (@var{rr}*@var{rr}) of doubles, covariance matrix of the structural innovations (noise in the state equation).
%! @item H
%! Vector (@var{pp}) of doubles, diagonal of covariance matrix of the measurement errors (corelation among measurement errors is handled by a model transformation).
%! @item Z
%! Matrix (@var{pp}*@var{mm}) of doubles, matrix relating the states to the observed variables.
%! @item mm
%! Integer scalar, number of state variables.
%! @item pp
%! Integer scalar, number of observed variables.
%! @item rr
%! Integer scalar, number of structural innovations.
%! @item Zflag
%! Integer scalar, equal to 0 if Z is a vector of indices targeting the obseved variables in the state vector, equal to 1 if Z is a @var{pp}*@var{mm} matrix.
%! @item diffuse_periods
%! Integer scalar, number of diffuse filter periods in the initialization step.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item dLIK
%! Double scalar, value of (minus) the likelihood for the first d observations.
%! @item likk
%! Column vector (d*1) of doubles, values of the density of each (fist) observations.
%! @item a
%! Vector (@var{mm}*1) of doubles, mean of the state vector at the end of the (sub)sample.
%! @item Pstar
%! Matrix (@var{mm}*@var{mm}) of doubles, covariance of the state vector at the end of the subsample.
%! @item llik
%! Matrix of doubles (d*@var{pp}) used by @ref{DsgeLikelihood_hh}
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{dsge_likelihood}, @ref{DsgeLikelihood_hh}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{univariate_kalman_filter_ss}
%! @end deftypefn
%@eod:
%
% Algorithm:
%
%   Uses the diffuse univariate filter as described in Durbin/Koopman (2012): "Time
%   Series Analysis by State Space Methods", Oxford University Press,
%   Second Edition, Ch. 5, 6.4 + 7.2.5

% Copyright © 2004-2024 Dynare Team
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

% Get sample size.
smpl = last-start+1;

% Initialize some variables.
dF   = 1;
isqvec = false;
if ndim(Q)>2
    Qvec = Q;
    Q=Q(:,:,1);
    isqvec = true;
end
QQ   = R*Q*transpose(R);   % Variance of R times the vector of structural innovations.
t    = start;              % Initialization of the time index.
dlikk= zeros(smpl,1);      % Initialization of the vector gathering the densities.
dLIK = Inf;                % Default value of the log likelihood.
oldK = Inf;
llik = zeros(smpl,pp);

newRank = rank(Pinf,diffuse_kalman_tol);
l2pi = log(2*pi);
s=0;

while newRank && (t<=last)
    s = t-start+1;
    d_index = data_index{t};
    for i=1:length(d_index)
        Zi = Z(d_index(i),:);
        prediction_error = Y(d_index(i),t) - Zi*a;      % nu_{t,i} in 6.13 in DK (2012)
        Fstar = Zi*Pstar*Zi' + H(d_index(i));           % F_{*,t} in 5.7 in DK (2012), relies on H being diagonal
        Finf  = Zi*Pinf*Zi';                            % F_{\infty,t} in 5.7 in DK (2012), relies on H being diagonal
        Kstar = Pstar*Zi';
        % Conduct check of rank
        % Pinf and Finf are always scaled such that their norm=1: Fstar/Pstar, instead,
        % depends on the actual values of std errors in the model and can be badly scaled.
        % experience is that diffuse_kalman_tol has to be bigger than kalman_tol, to ensure
        % exiting the diffuse filter properly, avoiding tests that provide false non-zero rank for Pinf.
        % Also the test for singularity is better set coarser for Finf than for Fstar for the same reason
        if Finf>diffuse_kalman_tol && newRank           % F_{\infty,t,i} = 0, use upper part of bracket on p. 175 DK (2012) for w_{t,i}
            Kinf   = Pinf*Zi';
            Kinf_Finf = Kinf/Finf;
            a         = a + Kinf_Finf*prediction_error;
            Pstar     = Pstar + Kinf*(Kinf_Finf'*(Fstar/Finf)) - Kstar*Kinf_Finf' - Kinf_Finf*Kstar';
            Pinf      = Pinf - Kinf*Kinf_Finf';
            llik(s,d_index(i)) = log(Finf) + l2pi;
            dlikk(s) = dlikk(s) + llik(s,d_index(i));
        elseif Fstar>kalman_tol
            llik(s,d_index(i)) = log(Fstar) + (prediction_error*prediction_error/Fstar) + l2pi;
            dlikk(s) = dlikk(s) + llik(s,d_index(i));
            a = a+Kstar*(prediction_error/Fstar);
            Pstar = Pstar-Kstar*(Kstar'/Fstar);
        else
            if Fstar<0 || Finf<0
                %pathological numerical case where variance is negative
                dLIK = NaN;
                return
            else
                % do nothing as a_{t,i+1}=a_{t,i} and P_{t,i+1}=P_{t,i}, see
                % p. 157, DK (2012)
            end
        end
    end
    if newRank
        oldRank = rank(Z*Pinf*Z',diffuse_kalman_tol);
    else
        oldRank = 0;
    end
    a     = T*a;
    if isqvec
        QQ = R*Qvec(:,:,t+1)*transpose(R);
    end
    Pstar = T*Pstar*T'+QQ;
    Pinf  = T*Pinf*T';
    if newRank
        newRank = rank(Z*Pinf*Z',diffuse_kalman_tol);
    end
    if oldRank ~= newRank
        disp('univariate_diffuse_kalman_filter:: T does influence the rank of Pinf!')
        disp('This may happen for models with order of integration >1.')
    end
    t = t+1;
end

if (t>last)
    warning('univariate_diffuse_kalman_filter:: There isn''t enough information to estimate the initial conditions of the nonstationary variables');
    dLIK = NaN;
    return
end

% Divide by two.
dlikk = .5*dlikk(1:s);
llik  = .5*llik(1:s,:);

dLIK = sum(dlikk(1+presample:end));
dlikk = llik;