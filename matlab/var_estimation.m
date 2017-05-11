function oo_ = var_estimation(M_, options_, oo_)
% function oo_ = var_estimation(M_, options_, oo_)
%
% INPUTS
%
%    M_:          [structure]  Model
%    options_:    [structure]  Options
%    oo_:         [structure]  Results
%
% OUTPUTS
%
%    oo_:         [structure]  Results

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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

%  Notation follows Lütkepohl (2005) Chapter 5

model = options_.var_estimation.model_name;
assert(isfield(M_.var, model));

if isfield(options_.var_estimation, 'datafile')
    datafile = options_.var_estimation.datafile;
else
    datafile = [model '.mat'];
end

var_list = cell(size(M_.var.(model).var_list_, 1), 1);
varMap = struct();
for i=1:size(M_.var.(model).var_list_, 1)
    var_list{i} = strtrim(M_.var.(model).var_list_(i,:));
    varMap.(var_list{i}) = i;
end
%varMap = containers.Map(var_list, 1:size(var_list,1));

p = M_.var.(model).order;
K = size(M_.var.(model).var_list_, 1);

%% Create Y and Z matrices from datafile
data = load(datafile);
T = data.estimationdata.nobs - p;
Y = zeros(K, T);
Z = ones(K*p+1, T);
for i = 1:K
    Y(i, :) = data.estimationdata.(var_list{i}).data(p+1:p+T)';
end

for i = 1:T
    Z(2:end, i) = vec(data.estimationdata{var_list{:}}.data(i+p-1:-1:i, :)');
end

%% Estimation
if ~isfield(M_.var.(model), 'restrictions')
    %% No restrictions
    B = Y*Z'/(Z*Z');
else
    %% Restrictions: Create C1, C2, c, R, r. Find gamma, B

    % Keep count of total restrictions in preprocessor
    N = 0; % The number of constraints
    for i = 1:size(M_.var.my_var_est.restrictions.exclusion_restrictions, 2)
        N = N + size(char(M_.var.my_var_est.restrictions.exclusion_restrictions{i}.restrictions{:,2}),1);
    end

    % Exclusion Restrictions
    C = zeros(N, K^2*p+K);
    c = zeros(N, 1);
    idx = 1;
    for i = 1:size(M_.var.my_var_est.restrictions.exclusion_restrictions, 2)
        lag = M_.var.my_var_est.restrictions.exclusion_restrictions{i}.lag;
        for j = 1:size(M_.var.my_var_est.restrictions.exclusion_restrictions{i}.restrictions, 2)
            eq = varMap.(M_.var.my_var_est.restrictions.exclusion_restrictions{i}.restrictions{j, 1});
            for k = 1:size(M_.var.my_var_est.restrictions.exclusion_restrictions{i}.restrictions{j, 2}, 1)
                varidx = varMap.(M_.var.my_var_est.restrictions.exclusion_restrictions{i}.restrictions{j, 2}(k,:));
                C(idx, (lag-1)*K^2+K + (varidx-1)*K + eq) = 1;
                idx = idx + 1;
            end
        end
    end

    % Variance Covariance Matrix restrictions
    %    if isfield(M_.var.(model).restrictions, 'covariance_const_restriction') ...
    %            || isfield(M_.var.(model).restrictions, 'covariance_pair_restriction')
    %    end

    % Find linear independent Columns of C
    [junk, beta1] = rref(C);
    beta = 1:K^2*p+K;
    beta1 = beta1';
    beta2 = setdiff(beta, beta1)';
    C1 = C(:, beta1);
    C2 = C(:, beta2);

    R = [C1\C2; eye(K^2*p+K-N)];
    r = C1\c;

    % Estimation...

end

oo_.var_estimation.(model).nu = B(:,1);
for i=1:K:K*p
    oo_.var_estimation.(model).A{(i-1)/K+1} = B(:, i+1:i+K);
end
end




