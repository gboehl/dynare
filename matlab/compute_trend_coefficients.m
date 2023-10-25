function [trend_addition, trend_coeff]=compute_trend_coefficients(M_,options_,nvarobs,ntobs)
% [trend_addition, trend_coeff]=compute_trend_coefficients(M_,options_,nvarobs,ntobs)
% Computes the trend coefficiencts and the trend, accounting for
% prefiltering
%
% INPUTS
%   M_              [structure] describing the model; called in the eval
%                               statement
%   options_        [structure] describing the options
%   nvarobs         [scalar]    number of observed variables
%   ntobs           [scalar]    length of data sample for estimation
%
% OUTPUTS
%   trend_addition  [nvarobs by ntobs double] matrix storing deterministic
%                               trend component
%   trend_coeff     [nvarobs by 1] vector storing trend slope
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2014-2023 Dynare Team
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


trend_coeff = zeros(nvarobs,1);
t = options_.trend_coeffs;
for i=1:length(t)
    if ~isempty(t{i})
        trend_coeff(i) = eval(t{i});
    end
end
trend_addition=trend_coeff*[options_.first_obs:options_.first_obs+ntobs-1];
if options_.prefilter
    trend_addition = bsxfun(@minus,trend_addition,mean(trend_addition,2));
end
