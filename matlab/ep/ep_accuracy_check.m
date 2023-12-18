function e = ep_accuracy_check(M_,options_,oo_)
% e = ep_accuracy_check(M_,options_,oo_)

% Copyright © 2016-2023 Dynare Team
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

endo_simul = oo_.endo_simul;
n = size(endo_simul,2);
[~, innovations, pfm, ~, ~, options_, oo_] = ...
    extended_path_initialization([], n-1, [], options_, M_, oo_);

options_.ep.accuracy.stochastic.order = options_.ep.stochastic.order;
[nodes,weights] = setup_integration_nodes(options_.ep.accuracy,pfm);

e = zeros(M_.endo_nbr,n);
for i=1:n
    e(:,i) = euler_equation_error(endo_simul(:,i),oo_.exo_simul, ...
                                  innovations, M_, options_,oo_,pfm,nodes,weights);
end
