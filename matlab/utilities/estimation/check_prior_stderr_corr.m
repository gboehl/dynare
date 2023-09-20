function check_prior_stderr_corr(estim_params_,bayestopt_)
% function check_prior_stderr_corr(estim_params_,bayestopt_)
% -------------------------------------------------------------------------
% Check if the prior allows for negative standard deviations and
% correlations larger than +-1. If so, issue a warning.
% -------------------------------------------------------------------------
% INPUTS
%  o estim_params_:           [struct] information on estimated parameters
%  o bayestopt_:              [struct] information on priors
% -------------------------------------------------------------------------
% OUTPUTS
%  none, but issues a warning if the prior allows for negative standard
% -------------------------------------------------------------------------
% This function is called by
%  o initial_estimation_checks.m
% -------------------------------------------------------------------------

% Copyright © 2023 Dynare Team
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

nvx = estim_params_.nvx; % number of stderr parameters for structural shocks
if nvx && any(bayestopt_.p3(1:nvx)<0)
    warning('Your prior allows for negative standard deviations for structural shocks. Due to working with variances, Dynare will be able to continue, but it is recommended to change your prior.')
end
offset = nvx;
nvn = estim_params_.nvn; % number of stderr parameters for measurement errors
if nvn && any(bayestopt_.p3(1+offset:offset+nvn)<0)
    warning('Your prior allows for negative standard deviations for measurement error. Due to working with variances, Dynare will be able to continue, but it is recommended to change your prior.')
end
offset = nvx+nvn;
ncx = estim_params_.ncx; % number of corr parameters for structural shocks
if ncx && (any(bayestopt_.p3(1+offset:offset+ncx)<-1) || any(bayestopt_.p4(1+offset:offset+ncx)>1))
    warning('Your prior allows for correlations between structural shocks larger than +-1 and will not integrate to 1 due to truncation. Please change your prior')
end
offset = nvx+nvn+ncx;
ncn = estim_params_.ncn; % number of corr parameters for measurement errors
if ncn && (any(bayestopt_.p3(1+offset:offset+ncn)<-1) || any(bayestopt_.p4(1+offset:offset+ncn)>1))
    warning('Your prior allows for correlations between measurement errors larger than +-1 and will not integrate to 1 due to truncation. Please change your prior')
end