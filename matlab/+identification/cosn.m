function [co, b, yhat] = cosn(H)
% function [co, b, yhat] = cosn(H)
% -------------------------------------------------------------------------
% computes the cosine of the angle between the (endogenous variable) H(:,1)
% and its projection onto the span of (exogenous variables) H(:,2:end)
% Note: This is not the same as multiple correlation coefficient since the
% means are not zero
% =========================================================================
% INPUTS
%   * H     [n by k]
%           Data matrix, endogenous variable y is in the first column,
%           exogenous variables X are in the remaining (k-1) columns
% -------------------------------------------------------------------------
% OUTPUTS
%   * co    [double] (approximate) multiple correlation coefficient
%   * b     [k by 1] ols estimator
%   * y     [n by 1] predicted endogenous values given ols estimation
% -------------------------------------------------------------------------
% This function is called by
%   * identification.checks.m
%   * ident_bruteforce.m
% =========================================================================
% Copyright © 2008-2019 Dynare Team
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
% =========================================================================

y = H(:,1);
X = H(:,2:end);

b=(X\y); %ols estimator
if any(isnan(b)) || any(isinf(b))
    b=0;
end
yhat =  X*b; %predicted values
if rank(yhat)
    co = abs(y'*yhat/sqrt((y'*y)*(yhat'*yhat)));
else
    co=0;
end
