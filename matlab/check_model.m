function check_model(M_)
% check_model(M_)
% Performs various consistency checks on the model

% Copyright (C) 2005-2023 Dynare Team
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

if ~all(M_.lead_lag_incidence(M_.maximum_lag+1,:) > 0)
    warning('Problem in the model specification: some variables don''t appear as current. Check whether this is desired.');
end

if any(diag(M_.Sigma_e)<0)
    error('You have specified negative shock variances. That is not allowed.')
end
if ~check_consistency_covariances(M_.Sigma_e)
    error('The specified covariances for the structural errors are not consistent with the variances as they imply a correlation larger than +-1')
end
if any(diag(M_.H)<0)
    error('You have specified negative measurement error variances. That is not allowed.')
end
if ~isequal(M_.H,0)
    if ~check_consistency_covariances(M_.H)
        error('The specified covariances for the measurement errors are not consistent with the variances as they imply a correlation larger than +-1')
    end
end