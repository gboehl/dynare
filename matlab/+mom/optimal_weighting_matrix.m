function W_opt = optimal_weighting_matrix(m_data, moments, q_lag)
% W_opt = optimal_weighting_matrix(m_data, moments, q_lag)
% -------------------------------------------------------------------------
% This function computes the optimal weigthing matrix by a Bartlett kernel with maximum lag q_lag
% Adapted from replication codes of Andreasen, Fernández-Villaverde, Rubio-Ramírez (2018):
% "The Pruned State-Space System for Non-Linear DSGE Models: Theory and Empirical Applications",
% Review of Economic Studies, 85(1):1-49.
% -------------------------------------------------------------------------
% INPUTS
%  o m_data                  [T x numMom]       selected data moments at each point in time
%  o moments                 [numMom x 1]       selected estimated moments (either data_moments or estimated model_moments)
%  o q_lag                   [integer]          Bartlett kernel maximum lag order
% -------------------------------------------------------------------------
% OUTPUTS 
%   o W_opt                  [numMom x numMom]  optimal weighting matrix
% -------------------------------------------------------------------------
% This function is called by
%  o mom.run.m
% -------------------------------------------------------------------------
% This function calls:
%  o corr_matrix (embedded)
% -------------------------------------------------------------------------

% Copyright © 2020-2023 Dynare Team
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


% initialize
[T,num_Mom] = size(m_data); % note that in m_data NaN values (due to leads or lags in matched_moments and missing data) were replaced by the mean

% center around moments (could be either data_moments or model_moments)
h_func = m_data - repmat(moments',T,1);

% the required correlation matrices
gamma_array = zeros(num_Mom,num_Mom,q_lag);
gamma0 = corr_matrix(h_func,T,num_Mom,0);
if q_lag > 0
    for ii=1:q_lag
        gamma_array(:,:,ii) = corr_matrix(h_func,T,num_Mom,ii);
    end
end

% the estimate of S
S = gamma0;
if q_lag > 0
    for ii=1:q_lag
        S = S + (1-ii/(q_lag+1))*(gamma_array(:,:,ii) + gamma_array(:,:,ii)');
    end
end

% the estimate of W
W_opt = S\eye(size(S,1));

W_opt = (W_opt+W_opt')/2; % ensure symmetry
end % main function end

% The correlation matrix
function gamma_corr = corr_matrix(h_func,T,num_Mom,v)
    gamma_corr = zeros(num_Mom,num_Mom);
    for t = 1+v:T
        gamma_corr = gamma_corr + h_func(t-v,:)'*h_func(t,:);
    end
    gamma_corr = gamma_corr/T;
end % corr_matrix end