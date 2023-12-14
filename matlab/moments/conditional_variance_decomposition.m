function [ConditionalVarianceDecomposition, ConditionalVarianceDecomposition_ME]= conditional_variance_decomposition(M_,options_,dr, Steps, SubsetOfVariables)
% [ConditionalVarianceDecomposition, ConditionalVarianceDecomposition_ME]= conditional_variance_decomposition(M_,options_,dr, Steps, SubsetOfVariables)
% This function computes the conditional variance decomposition of a given state space model
% for a subset of endogenous variables.
%
% INPUTS
%   M_                  [struct]        Matlab's structure describing the Model (initialized by dynare, see @ref{M_}).          
%   options_            [struct]        Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%   dr :                [struct]        Dynare decision rules structure
%   Steps               [integer]       1*h vector of dates.
%   SubsetOfVariables   [integer]       1*q vector of indices (declaration order).
%
% OUTPUTS
%   ConditionalVarianceDecomposition  [double] [n h p] array, where
%                                                    n is equal to length(SubsetOfVariables)
%                                                    h is the number of Steps
%                                                    p is the number of state innovations and
%   ConditionalVarianceDecomposition_ME  [double] [m h p] array, where
%                                                    m is equal to length(intersect(SubsetOfVariables,varobs))
%                                                    h is the number of Steps
%                                                    p is the number of state innovations and

% Copyright © 2010-2023 Dynare Team
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

if any(Steps <= 0)
    error(['Conditional variance decomposition: All periods must be strictly ' ...
           'positive'])
end

[transition_matrix, impulse_matrix] = kalman_transition_matrix(dr,(1:M_.endo_nbr)',M_.nstatic+(1:M_.nspred)');

number_of_state_innovations = M_.exo_nbr;
number_of_state_equations = M_.endo_nbr;
order_var = dr.order_var;
nSteps = length(Steps);

ConditionalVariance = zeros(number_of_state_equations,nSteps,number_of_state_innovations);

if M_.sigma_e_is_diagonal
    B = impulse_matrix.* repmat(sqrt(diag(M_.Sigma_e)'),number_of_state_equations,1);
else
    B = impulse_matrix*chol(M_.Sigma_e)';
end

for i=1:number_of_state_innovations
    BB = B(:,i)*B(:,i)';
    V = zeros(number_of_state_equations,number_of_state_equations);
    m = 1;
    for h = 1:max(Steps)
        V = transition_matrix*V*transition_matrix'+BB;
        if h == Steps(m)
            ConditionalVariance(order_var,m,i) = diag(V);
            m = m+1;
        end
    end
end

ConditionalVariance = ConditionalVariance(SubsetOfVariables,:,:);

NumberOfVariables = length(SubsetOfVariables);
SumOfVariances = zeros(NumberOfVariables,nSteps);
for h = 1:length(Steps)
    SumOfVariances(:,h) = sum(ConditionalVariance(:,h,:),3);
end

ConditionalVarianceDecomposition = zeros(NumberOfVariables,length(Steps),number_of_state_innovations);
for i=1:number_of_state_innovations
    for h = 1:length(Steps)
        ConditionalVarianceDecomposition(:,h,i) = squeeze(ConditionalVariance(:,h,i))./SumOfVariances(:,h);
    end
end

% get intersection of requested variables and observed variables with
% Measurement error

if ~all(diag(M_.H)==0)
    if isoctave && octave_ver_less_than('8.4') %Octave bug #60347
        [observable_pos,index_subset,index_observables]=intersect_stable(SubsetOfVariables,options_.varobs_id);
    else
        [observable_pos,index_subset,index_observables]=intersect(SubsetOfVariables,options_.varobs_id,'stable');
    end
    ME_Variance=diag(M_.H);

    ConditionalVarianceDecomposition_ME = zeros(length(observable_pos),length(Steps),number_of_state_innovations+1);
    for i=1:number_of_state_innovations
        for h = 1:length(Steps)
            ConditionalVarianceDecomposition_ME(:,h,i) = squeeze(ConditionalVariance(index_subset,h,i))./(SumOfVariances(index_subset,h)+ME_Variance(index_observables));
        end
    end
    ConditionalVarianceDecomposition_ME(:,:,number_of_state_innovations+1)=1-sum(ConditionalVarianceDecomposition_ME(:,:,1:number_of_state_innovations),3);
else
    ConditionalVarianceDecomposition_ME=[];
end
