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
p = length(oo_.var.(var_model_name).ar);

% Get the number of variables
n = length(oo_.var.(var_model_name).ar{1});

% Initialise the companion matrix
oo_.var.(var_model_name).H = zeros(n*p);

% Fill the companion matrix
oo_.var.(var_model_name).H(1:n,1:n) = oo_.var.(var_model_name).ar{1};

if p>1
    for i=2:p
        oo_.var.(var_model_name).H(1:n,(i-1)*n+(1:n)) = oo_.var.(var_model_name).ar{i};
        oo_.var.(var_model_name).H((i-1)*n+(1:n),(i-2)*n+(1:n)) = eye(n); 
    end
end