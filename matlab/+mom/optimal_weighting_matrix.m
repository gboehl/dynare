function W_opt = optimal_weighting_matrix(m_data, moments, q_lag)
% W_opt = optimal_weighting_matrix(m_data, moments, q_lag)
% -------------------------------------------------------------------------
% This function computes the optimal weigthing matrix by a Bartlett kernel with maximum lag q_lag
% Adapted from replication codes of
%  o Andreasen, Fernández-Villaverde, Rubio-Ramírez (2018): "The Pruned State-Space System for Non-Linear DSGE Models: Theory and Empirical Applications", Review of Economic Studies, 85(1):1-49.
% =========================================================================
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
%  o CorrMatrix (embedded)
% =========================================================================
% Copyright (C) 2020-2021 Dynare Team
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
% -------------------------------------------------------------------------
% Author(s): 
% o Willi Mutschler (willi@mutschler.eu)
% o Johannes Pfeifer (jpfeifer@uni-koeln.de)
% =========================================================================

% Initialize
[T,num_Mom] = size(m_data); %note that in m_data NaN values (due to leads or lags in matched_moments and missing data) were replaced by the mean

% center around moments (could be either data_moments or model_moments)
h_Func = m_data - repmat(moments',T,1);

% The required correlation matrices
GAMA_array = zeros(num_Mom,num_Mom,q_lag);
GAMA0 = Corr_Matrix(h_Func,T,num_Mom,0);
if q_lag > 0
    for ii=1:q_lag
        GAMA_array(:,:,ii) = Corr_Matrix(h_Func,T,num_Mom,ii);
    end
end

% The estimate of S
S = GAMA0;
if q_lag > 0
    for ii=1:q_lag
        S = S + (1-ii/(q_lag+1))*(GAMA_array(:,:,ii) + GAMA_array(:,:,ii)');
    end
end

% The estimate of W
W_opt = S\eye(size(S,1));

end

% The correlation matrix
function GAMA_corr = Corr_Matrix(h_Func,T,num_Mom,v)
    GAMA_corr = zeros(num_Mom,num_Mom);
    for t = 1+v:T
        GAMA_corr = GAMA_corr + h_Func(t-v,:)'*h_Func(t,:);
    end
    GAMA_corr = GAMA_corr/T;
end
