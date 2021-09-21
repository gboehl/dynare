function oo_=fill_mh_mode(xparam1,stdh,M_,options_,estim_params_,bayestopt_,oo_, field_name)
%function oo_=fill_mh_mode(xparam1,stdh,M_,options_,estim_params_,bayestopt_,oo_, field_name)
%
% INPUTS
%   o xparam1       [double]   (p*1) vector of estimate parameters.
%   o stdh          [double]   (p*1) vector of estimate parameters.
%   o M_                        Matlab's structure describing the Model (initialized by dynare, see @ref{M_}).
%   o estim_params_             Matlab's structure describing the estimated_parameters (initialized by dynare, see @ref{estim_params_}).
%   o options_                  Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%   o bayestopt_                Matlab's structure describing the priors (initialized by dynare, see @ref{bayesopt_}).
%   o oo_                       Matlab's structure gathering the results (initialized by dynare, see @ref{oo_}).
%
% OUTPUTS
%   o oo_                       Matlab's structure gathering the results
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 2005-2021 Dynare Team
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

nvx = estim_params_.nvx;  % Variance of the structural innovations (number of parameters).
nvn = estim_params_.nvn;  % Variance of the measurement innovations (number of parameters).
ncx = estim_params_.ncx;  % Covariance of the structural innovations (number of parameters).
ncn = estim_params_.ncn;  % Covariance of the measurement innovations (number of parameters).
np  = estim_params_.np ;  % Number of deep parameters.

if np
    ip = nvx+nvn+ncx+ncn+1;
    for i=1:np
        name = bayestopt_.name{ip};
        oo_.([field_name '_mode']).parameters.(name) = xparam1(ip);
        oo_.([field_name '_std_at_mode']).parameters.(name) = stdh(ip);
        ip = ip+1;
    end
end
if nvx
    ip = 1;
    for i=1:nvx
        k = estim_params_.var_exo(i,1);
        name = M_.exo_names{k};
        oo_.([field_name '_mode']).shocks_std.(name)= xparam1(ip);
        oo_.([field_name '_std_at_mode']).shocks_std.(name) = stdh(ip);
        ip = ip+1;
    end
end
if nvn
    ip = nvx+1;
    for i=1:nvn
        name = options_.varobs{estim_params_.nvn_observable_correspondence(i,1)};
        oo_.([field_name '_mode']).measurement_errors_std.(name) = xparam1(ip);
        oo_.([field_name '_std_at_mode']).measurement_errors_std.(name) = stdh(ip);
        ip = ip+1;
    end
end

if ncx
    ip = nvx+nvn+1;
    for i=1:ncx
        k1 = estim_params_.corrx(i,1);
        k2 = estim_params_.corrx(i,2);
        NAME = [M_.exo_names{k1} '_' M_.exo_names{k2}];
        oo_.([field_name '_mode']).shocks_corr.(name) = xparam1(ip);
        oo_.([field_name '_std_at_mode']).shocks_corr.(name) = stdh(ip);
        ip = ip+1;
    end
end

if ncn
    ip = nvx+nvn+ncx+1;
    for i=1:ncn
        k1 = estim_params_.corrn(i,1);
        k2 = estim_params_.corrn(i,2);
        NAME = [M_.endo_names{k1} '_' M_.endo_names{k2}];
        oo_.([field_name '_mode']).measurement_errors_corr.(name) = xparam1(ip);
        oo_.([field_name '_std_at_mode']).measurement_errors_corr.(name) = stdh(ip);
        ip = ip+1;
    end
end