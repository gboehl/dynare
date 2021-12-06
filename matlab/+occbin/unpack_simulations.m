function oo_=unpack_simulations(M_,oo_,options_)
% function oo_=unpack_simulations(M_,oo_,options_)
% Writes Occbin simulations from matrix to structure
% 
% Inputs
% - M_                  [structure]     Matlab's structure describing the model
% - oo_                 [structure]     Matlab's structure containing the results
% - options_            [structure]     Matlab's structure containing the options
%
% Outputs
% - oo_                 [structure]     Matlab's structure containing the results

% Copyright (C) 2021 Dynare Team
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

for i=1:M_.endo_nbr
    % unpack the IRFs
    oo_.occbin.endo_linear.(M_.endo_names{i})= oo_.occbin.simul.linear(:,i);
    oo_.occbin.endo_piecewise.(M_.endo_names{i})=oo_.occbin.simul.piecewise(:,i);
    oo_.occbin.endo_ss.(M_.endo_names{i})=oo_.occbin.simul.ys(i);
end
for i=1:length(oo_.occbin.simul.exo_pos)
    oo_.occbin.exo.(M_.exo_names{i})=options_.occbin.simul.SHOCKS(:,i);
end