function estim_params_= get_matrix_entries_for_psd_check(M_,estim_params_)
% function estim_params_= get_matrix_entries_for_psd_check(M_)
% Get entries of Sigma_e and H to check for positive definiteness
%
% INPUTS
%   M_:             structure storing the model information
%   estim_params_:  structure storing information about estimated
%                   parameters
% OUTPUTS
%   estim_params_:  structure storing information about estimated
%                   parameters
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2020 Dynare Team
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

%% get the non-zero rows and columns of Sigma_e and H

H_non_zero_rows=find(~all(M_.H==0,1));
H_non_zero_columns=find(~all(M_.H==0,2));
if ~isequal(H_non_zero_rows,H_non_zero_columns') || (any(any(M_.H-M_.H'>1e-10)))
    error('Measurement error matrix not symmetric')
end
if isfield(estim_params_,'nvn_observable_correspondence')
    estim_params_.H_entries_to_check_for_positive_definiteness=union(H_non_zero_rows,estim_params_.nvn_observable_correspondence(:,1));
else
    estim_params_.H_entries_to_check_for_positive_definiteness=H_non_zero_rows;
end

Sigma_e_non_zero_rows=find(~all(M_.Sigma_e==0,1));
Sigma_e_non_zero_columns=find(~all(M_.Sigma_e==0,2));
if ~isequal(Sigma_e_non_zero_rows,Sigma_e_non_zero_columns') || (any(any(M_.Sigma_e-M_.Sigma_e'>1e-10)))
    error('Structual error matrix not symmetric')
end
if isfield(estim_params_,'var_exo') && ~isempty(estim_params_.var_exo)
    estim_params_.Sigma_e_entries_to_check_for_positive_definiteness=union(Sigma_e_non_zero_rows,estim_params_.var_exo(:,1));
else
    estim_params_.Sigma_e_entries_to_check_for_positive_definiteness=Sigma_e_non_zero_rows;
end