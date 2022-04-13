function moments=compute_model_moments(dr,M_,options_)
%
% INPUTS
%    dr:        structure describing model solution
%    M_:   structure of Dynare options
%     options_
%
% OUTPUTS
%    moments: a cell array containing requested moments
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2008-2017 Dynare Team
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

[ivar,vartan,options_] = get_variables_list(options_,M_);
moments = th_autocovariances(dr,ivar,M_,options_,options_.nodecomposition);
