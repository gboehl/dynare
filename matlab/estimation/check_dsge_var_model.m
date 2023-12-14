function check_dsge_var_model(M_, estim_params_, bayestopt_)

% Check if the dsge model can be estimated with the DSGE-VAR approach.

% Copyright © 2013-2014 Dynare Team
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

if estim_params_.nvn
    error('Estimation::DsgeVarLikelihood: Measurement errors are not allowed!')
end

if estim_params_.ncn
    error('Estimation::DsgeVarLikelihood: Measurement errors are not allowed!')
end

if any(vec(M_.H))
    error('Estimation::DsgeVarLikelihood: Measurement errors are not allowed!')
end

if estim_params_.ncx
    error('Estimation::DsgeVarLikelihood: Structural innovations cannot be correlated using Dynare''s interface! Introduce the correlations in the model block instead.')
end

if M_.exo_nbr>1 && any(vec(tril(M_.Sigma_e,-1)))
    error('Estimation::DsgeVarLikelihood: Structural innovations cannot be correlated using Dynare''s interface! Introduce the correlations in the model block instead.')
end

if isequal(bayestopt_.with_trend,1)
    error('Estimation::DsgeVarLikelihood: Linear trend is not yet implemented!')
end
