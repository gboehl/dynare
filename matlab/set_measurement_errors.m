function M_ = set_measurement_errors(xparam1,estim_params_,M_)
% function M_=set_measurement_errors(xparam1,estim_params_,M_)
% Sets parameters value (except measurement errors)
% This is called for computations such as IRF and forecast
% when measurement errors aren't taken into account; in contrast to
% set_parameters.m, the global M_-structure is not altered
%
% INPUTS
%    xparam1:   vector of parameters to be estimated (initial values)
%    M_:        Dynare model-structure
%
% OUTPUTS
%    M_:        Dynare model-structure
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2017 Dynare Team
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

H = M_.H;
Correlation_matrix_ME = M_.Correlation_matrix_ME;

% setting measument error variance; on the diagonal of Covariance matrix; used later
% for updating covariances
offset = estim_params_.nvx;

if estim_params_.nvn
    for i=1:estim_params_.nvn
        k = estim_params_.nvn_observable_correspondence(i,1);
        H(k,k) = xparam1(i+offset)^2;
    end
end

% update offset
offset = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx;

% setting measurement error covariances
if estim_params_.ncn
    corrn_observable_correspondence = estim_params_.corrn_observable_correspondence;
    for i=1:estim_params_.ncn
        k1 = corrn_observable_correspondence(i,1);
        k2 = corrn_observable_correspondence(i,2);
        Correlation_matrix_ME(k1,k2) = xparam1(i+offset);
        Correlation_matrix_ME(k2,k1) = Correlation_matrix_ME(k1,k2);
    end
end
%build covariance matrix from correlation matrix and variances already on
%diagonal
H = diag(sqrt(diag(H)))*Correlation_matrix_ME*diag(sqrt(diag(H)));
%if calibrated covariances, set them now to their stored value
if isfield(estim_params_,'calibrated_covariances_ME')
    H(estim_params_.calibrated_covariances_ME.position)=estim_params_.calibrated_covariances_ME.cov_value;
end

M_.H = H;
M_.Correlation_matrix_ME=Correlation_matrix_ME;
