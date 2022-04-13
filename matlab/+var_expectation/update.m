function update(varexpectationmodelname)

% Updates the parameters of a VAR_EXPECTATION_MODEL.
%
% INPUTS
% - varepxpectationmodelname       [string]    Name of the VAR expectation model.
% 
% OUTPUTS
% None

% Copyright © 2018-2021 Dynare Team
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

global M_ oo_

M_ = var_expectation.update_parameters(varexpectationmodelname, M_, oo_);
