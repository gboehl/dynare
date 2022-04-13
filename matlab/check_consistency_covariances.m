function correct_flag=check_consistency_covariances(Covariance_matrix)
% function check_consistency_covariances(Covariance_matrix)
% checks consistency of covariance matrices by checking whether the
% covariances imply correlations bigger than 1.
%
% Outputs: correct_flag         [scalar] 0 if not consistent, 1 otherwise
% Inputs: Covariance_matrix     [matrix] covariance matrix to be checked

% Copyright © 2013 Dynare Team
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

%compute theoretical bound by assuming correlation of 1
bound=diag(sqrt(diag(Covariance_matrix)))*ones(size(Covariance_matrix))*diag(sqrt(diag(Covariance_matrix)));
correct_flag=1;
% We must check that Covariance_matrix is non-empty, because tril([],-1) fails under Octave
if length(Covariance_matrix) > 0 && (any(any(tril(Covariance_matrix,-1)>bound)) || any(any(tril(Covariance_matrix,-1)<-bound)))
    correct_flag=0;
end
