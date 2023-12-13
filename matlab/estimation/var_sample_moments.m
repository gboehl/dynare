function dataset_info=var_sample_moments(nlag, var_trend_order, dataset_, dataset_info)
% dataset_info=var_sample_moments(nlag, var_trend_order, dataset_, dataset_info)
% Computes the sample moments of a VAR model.
%
% The VAR(p) model is defined by:
%
%   y_t = \sum_{k=1}^p y_{t-k} A_k + z_t C + e_t  for t = 1,...,T
%
% where y_t is a 1*m vector of observed endogenous variables, p is the
% number of lags, A_k is an m*m real matrix, z_t is a 1*q vector of
% exogenous (deterministic) variables, C is a q*m real matrix and
% e_t is a vector of exogenous stochastic shocks. T is the number
% of observations. The deterministic exogenous variables are assumed to
% be a polynomial trend of order q = "var_trend_order".
%
% We define:
%
%  <>  Y = (y_1',y_2',...,y_T')' a T*m matrix,
%
%  <>  x_t = (y_{t-1},y_{t-2},...,y_{t-p},z_t) a 1*(mp+q) row vector,
%
%  <>  X = (x_1',x_2',...,x_T')' a T*(mp+q) matrix,
%
%  <>  E = (e_1',e_2',...,e_T')' a T*m matrix and
%
%  <>  A = (A_1',A_2',...,A_p',C')' an (mp+q)*m matrix of coefficients.
%
% So that we can equivalently write the VAR(p) model using the following
% matrix representation:
%
%   Y = X * A +E
%
%
% INPUTS
%   o nlag                [integer] Number of lags in the VAR model.
%   o var_trend_order     [integer] Order of the polynomial exogenous trend:
%                                       = -1 no constant and no linear trend,
%                                       =  0 constant and no linear trend,
%                                       =  1 constant and linear trend.
%   o dataset_            [dseries] The sample.
%   o dataset_info        [structure] data set information
%
% OUTPUTS
%   dataset_info          [structure] containing the following new fields
%   o mYY                 [double]  Y'*Y an m*m matrix.
%   o mXY                 [double]  X'*Y an (mp+q)*m matrix.
%   o mYX                 [double]  Y'*X an m*(mp+q) matrix.
%   o mXX                 [double]  X'*X an (mp+q)*(mp+q) matrix.
%   o Ydata               [double]  Y a T*m matrix.
%   o Xdata               [double]  X a T*(mp+q) matrix.
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright © 2007-2023 Dynare Team
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

LastObservation = dataset_.dates(end);
FirstObservation = dataset_.dates(1)+nlag;

NumberOfObservations = LastObservation-FirstObservation+1;
NumberOfVariables = dataset_.vobs;

if isequal(var_trend_order,-1)
    % No constant no linear trend case.
    X = zeros(NumberOfObservations,NumberOfVariables*nlag);
elseif isequal(var_trend_order,0)
    % Constant and no linear trend case.
    X = ones(NumberOfObservations,NumberOfVariables*nlag+1);
elseif isequal(var_trend_order,1)
    % Constant and linear trend case.
    X = ones(NumberOfObservations,NumberOfVariables*nlag+2);
else
    error('Estimation::var_sample_moments: trend must be equal to -1,0 or 1!')
end

% I build matrices Y and X
Y = dataset_(FirstObservation:LastObservation).data;

for t=1:NumberOfObservations
    currentdate = FirstObservation+(t-1);
    for lag = 1:nlag
        X(t,(lag-1)*NumberOfVariables+1:lag*NumberOfVariables) = dataset_(currentdate-lag).data;
    end
end

if (var_trend_order == 1)
    X(:,end) = transpose(1:NumberOfObservations);
end
dataset_info.mYY=Y'*Y;
dataset_info.mYX=Y'*X;
dataset_info.mXY=X'*Y;
dataset_info.mXX=X'*X;
dataset_info.Ydata=Y;
dataset_info.Xdata=X;