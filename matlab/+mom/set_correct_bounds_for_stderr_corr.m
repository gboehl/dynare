function BoundsInfo = set_correct_bounds_for_stderr_corr(estim_params_,BoundsInfo)
% BoundsInfo = set_correct_bounds_for_stderr_corr(estim_params_,BoundsInfo)
% -------------------------------------------------------------------------
% Set correct bounds for standard deviation and corrrelation parameters
% -------------------------------------------------------------------------
% INPUTS
% o estim_params_ [struct] information on estimated parameters
% o BoundsInfo    [struct] information on bounds
% -------------------------------------------------------------------------
% OUTPUT
% o BoundsInfo    [struct] updated bounds
% -------------------------------------------------------------------------
% This function is called by
%  o mom.run
% -------------------------------------------------------------------------

% Copyright © 2023 Dynare Team
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


number_of_estimated_parameters = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.np;
% set correct bounds for standard deviations and corrrelations
param_of_interest = (1:number_of_estimated_parameters)'<=estim_params_.nvx+estim_params_.nvn;
LB_below_0 = (BoundsInfo.lb<0 & param_of_interest);
BoundsInfo.lb(LB_below_0) = 0;
param_of_interest = (1:number_of_estimated_parameters)'> estim_params_.nvx+estim_params_.nvn & (1:number_of_estimated_parameters)'<estim_params_.nvx+estim_params_.nvn +estim_params_.ncx + estim_params_.ncn;
LB_below_minus_1 = (BoundsInfo.lb<-1 & param_of_interest);
UB_above_1 = (BoundsInfo.ub>1 & param_of_interest);
BoundsInfo.lb(LB_below_minus_1) = -1;
BoundsInfo.ub(UB_above_1) = 1;