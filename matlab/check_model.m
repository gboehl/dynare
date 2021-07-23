function check_model(DynareModel)

% Copyright (C) 2005-2013 Dynare Team
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

if ~all(DynareModel.lead_lag_incidence(DynareModel.maximum_lag+1,:) > 0)
    warning('Problem in the model specification: some variables don''t appear as current. Check whether this is desired.');
end

if ~check_consistency_covariances(DynareModel.Sigma_e)
    error('The specified covariances for the structural errors are not consistent with the variances as they imply a correlation larger than +-1')
end
if ~isequal(DynareModel.H,0)
    if ~check_consistency_covariances(DynareModel.H)
        error('The specified covariances for the measurement errors are not consistent with the variances as they imply a correlation larger than +-1')
    end
end