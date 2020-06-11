function Wopt = method_of_moments_optimal_weighting_matrix(m_data, moments, qLag)
% Wopt = method_of_moments_optimal_weighting_matrix(m_data, moments, qLag)
% -------------------------------------------------------------------------
% This function computes the optimal weigthing matrix by a Bartlett kernel with maximum lag qlag
% Adapted from replication codes of
%  o Andreasen, Fernández-Villaverde, Rubio-Ramírez (2018): "The Pruned State-Space System for Non-Linear DSGE Models: Theory and Empirical Applications", Review of Economic Studies, 85(1):1-49.
% =========================================================================
% INPUTS
%  o m_data                  [T x numMom]       selected empirical or theoretical moments at each point in time
%  o moments                 [numMom x 1]       mean of selected empirical or theoretical moments
%  o qlag                    [integer]          Bartlett kernel maximum lag order
% -------------------------------------------------------------------------
% OUTPUTS 
%   o Wopt                   [numMom x numMom]  optimal weighting matrix
% -------------------------------------------------------------------------
% This function is called by
%  o method_of_moments.m
% -------------------------------------------------------------------------
% This function calls:
%  o CorrMatrix (embedded)
% =========================================================================
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------
% Author(s): 
% o Willi Mutschler (willi@mutschler.eu)
% o Johannes Pfeifer (jpfeifer@uni-koeln.de)
% =========================================================================

% Initialize
[T,numMom] = size(m_data); %note that in m_data nan values (due to leads or lags in matchedmoments) are removed so T is the effective sample size

% center around moments (could be either datamoments or modelmoments)
hFunc = m_data - repmat(moments',T,1);

% The required correlation matrices
GAMA_array = zeros(numMom,numMom,qLag);
GAMA0 = CorrMatrix(hFunc,T,numMom,0);
if qLag > 0
    for ii=1:qLag
        GAMA_array(:,:,ii) = CorrMatrix(hFunc,T,numMom,ii);
    end
end

% The estimate of S
S = GAMA0;
if qLag > 0
    for ii=1:qLag
        S = S + (1-ii/(qLag+1))*(GAMA_array(:,:,ii) + GAMA_array(:,:,ii)');
    end
end

% The estimate of W
Wopt = S\eye(size(S,1));

% Check positive definite W
try 
    chol(Wopt);
catch err
    error('method_of_moments: The optimal weighting matrix is not positive definite. Check whether your model implies stochastic singularity\n')
end

end

% The correlation matrix
function GAMAcorr = CorrMatrix(hFunc,T,numMom,v)
    GAMAcorr = zeros(numMom,numMom);
    for t = 1+v:T
        GAMAcorr = GAMAcorr + hFunc(t-v,:)'*hFunc(t,:);
    end
    GAMAcorr = GAMAcorr/T;
end