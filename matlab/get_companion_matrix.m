function get_companion_matrix(var_model_name)

% Gets the companion matrix associated with the var specified by
% var_model_name. Output stored in cellarray oo_.var.(var_model_name).H.
%
% INPUTS
% - var_model_name   [string]        the name of the VAR model
%
% OUTPUTS
% - None

% Copyright (C) 2018 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global oo_

get_ar_matrices(var_model_name);

% Get the number of lags
p = length(oo_.var.(var_model_name).AutoregressiveMatrices);

% Get the number of variables
n = length(oo_.var.(var_model_name).AutoregressiveMatrices{1});

if all(cellfun(@iszero, oo_.var.(var_model_name).ecm))
    % Build the companion matrix (standard VAR)
    oo_.var.(var_model_name).CompanionMatrix = zeros(n*p);
    oo_.var.(var_model_name).CompanionMatrix(1:n,1:n) = oo_.var.(var_model_name).AutoregressiveMatrices{1};
    if p>1
        for i=2:p
            oo_.var.(var_model_name).CompanionMatrix(1:n,(i-1)*n+(1:n)) = oo_.var.(var_model_name).AutoregressiveMatrices{i};
            oo_.var.(var_model_name).CompanionMatrix((i-1)*n+(1:n),(i-2)*n+(1:n)) = eye(n);
        end
    end
else
    B = zeros(n,n,p+1);
    idx = oo_.var.(var_model_name).ecm_idx;
    B(:,:,1) = oo_.var.(var_model_name).AutoregressiveMatrices{1};
    B(idx, idx, 1) = B(idx,idx, 1) + eye(length(idx));
    for i=2:p
        B(idx,idx,i) = oo_.var.(var_model_name).AutoregressiveMatrices{i}(idx,idx)-oo_.var.(var_model_name).AutoregressiveMatrices{i-1}(idx,idx);
    end
    B(idx,idx,p+1) = -oo_.var.(var_model_name).AutoregressiveMatrices{p}(idx,idx)
    % Build the companion matrix (VECM, rewrite in levels)
    oo_.var.(var_model_name).CompanionMatrix = zeros(n*(p+1));
    for i=1:p
        oo_.var.(var_model_name).CompanionMatrix(1:n, (i-1)*n+(1:n)) = B(:,:,1);
        oo_.var.(var_model_name).CompanionMatrix(i*n+(1:n),(i-1)*n+(1:n)) = eye(n);
    end
    oo_.var.(var_model_name).CompanionMatrix(1:n, p*n+(1:n)) = B(:,:,p+1);
end