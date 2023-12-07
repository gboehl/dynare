function M_ = set_all_parameters(xparam1,estim_params_,M_)

%@info:
%! @deftypefn {Function File} {@var{M_} =} dseries (@var{xparams1},@var{estim_params_},@var{M_})
%! @anchor{set_all_parameters}
%! @sp 1
%! Update parameter values (deep parameters and covariance matrices).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item xparam1
%! N*1 vector of doubles, the values of the N estimated parameters.
%! @item estim_params_
%! Dynare structure describing the estimated parameters.
%! @item M_
%! Dynare structure describing the model.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item M_
%! Dynare structure describing the model, with updated parameters and covariances matrices.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{DsgeSmoother}, @ref{dynare_estimation_1}, @ref{@@gsa.monte_carlo_filtering}, @ref{identification.analysis}, @ref{PosteriorFilterSmootherAndForecast}, @ref{prior_posterior_statistics_core}, @ref{prior_sampler}
%! @sp 2
%! @strong{This function calls:}
%! @sp 2
%! @end deftypefn
%@eod:

% Copyright © 2003-2017 Dynare Team
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

nvx = estim_params_.nvx;
ncx = estim_params_.ncx;
nvn = estim_params_.nvn;
ncn = estim_params_.ncn;
np = estim_params_.np;
if nvx || ncx
    Sigma_e = M_.Sigma_e;
    Correlation_matrix = M_.Correlation_matrix;
end
H = M_.H;
Correlation_matrix_ME = M_.Correlation_matrix_ME;
% setting shocks variance on the diagonal of Covariance matrix; used later
% for updating covariances
if nvx
    var_exo = estim_params_.var_exo;
    for i=1:nvx
        k =var_exo(i,1);
        Sigma_e(k,k) = xparam1(i)^2;
    end
end
% update offset
offset = nvx;

% setting measument error variance; on the diagonal of Covariance matrix; used later
% for updating covariances
if nvn
    for i=1:nvn
        k = estim_params_.nvn_observable_correspondence(i,1);
        H(k,k) = xparam1(i+offset)^2;
    end
end

% update offset
offset = nvx+nvn;

% setting shocks covariances
if ncx
    corrx = estim_params_.corrx;
    for i=1:ncx
        k1 = corrx(i,1);
        k2 = corrx(i,2);
        Correlation_matrix(k1,k2) = xparam1(i+offset);
        Correlation_matrix(k2,k1) = Correlation_matrix(k1,k2);
    end
end
% update offset
offset = nvx+nvn+ncx;

% setting measurement error covariances
if ncn
    corrn_observable_correspondence = estim_params_.corrn_observable_correspondence;
    for i=1:ncn
        k1 = corrn_observable_correspondence(i,1);
        k2 = corrn_observable_correspondence(i,2);
        Correlation_matrix_ME(k1,k2) = xparam1(i+offset);
        Correlation_matrix_ME(k2,k1) = Correlation_matrix_ME(k1,k2);
    end
end

% update offset
offset = nvx+ncx+nvn+ncn;
% setting structural parameters
%
if np
    M_.params(estim_params_.param_vals(:,1)) = xparam1(offset+1:end);
end

% updating matrices in M_
if nvx || ncx
    %build covariance matrix from correlation matrix and variances already on
    %diagonal
    Sigma_e = diag(sqrt(diag(Sigma_e)))*Correlation_matrix*diag(sqrt(diag(Sigma_e)));
    %if calibrated covariances, set them now to their stored value
    if isfield(estim_params_,'calibrated_covariances')
        Sigma_e(estim_params_.calibrated_covariances.position)=estim_params_.calibrated_covariances.cov_value;
    end
    M_.Sigma_e = Sigma_e;
    M_.Correlation_matrix=Correlation_matrix;
end
if nvn || ncn
    %build covariance matrix from correlation matrix and variances already on
    %diagonal
    H = diag(sqrt(diag(H)))*Correlation_matrix_ME*diag(sqrt(diag(H)));
    %if calibrated covariances, set them now to their stored value
    if isfield(estim_params_,'calibrated_covariances_ME')
        H(estim_params_.calibrated_covariances_ME.position)=estim_params_.calibrated_covariances_ME.cov_value;
    end
    M_.H = H;
    M_.Correlation_matrix_ME=Correlation_matrix_ME;
end