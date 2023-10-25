function [shocks, spfm_exo_simul, innovations, oo_] = extended_path_shocks(innovations, ep, exogenousvariables, sample_size,M_,options_, oo_)
% [shocks, spfm_exo_simul, innovations, oo_] = extended_path_shocks(innovations, ep, exogenousvariables, sample_size,M_,options_, oo_)
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

% Simulate shocks.
if isempty(exogenousvariables)
    switch ep.innovation_distribution
      case 'gaussian'
        shocks = zeros(sample_size, M_.exo_nbr);
        shocks(:,innovations.positive_var_indx) = transpose(transpose(innovations.covariance_matrix_upper_cholesky)*randn(innovations.effective_number_of_shocks,sample_size));
      case 'calibrated'
        options = options_;
        options.periods = options.ep.periods;
        oo_local = make_ex_(M_,options,oo_);
        shocks = oo_local.exo_simul(2:end,:);
      otherwise
        error(['extended_path:: ' ep.innovation_distribution ' distribution for the structural innovations is not (yet) implemented!'])
    end
else
    shocks = exogenousvariables;
    innovations.positive_var_indx = find(sum(abs(shocks)>0));
end

% Copy the shocks in exo_simul
oo_.exo_simul = shocks;
spfm_exo_simul = repmat(oo_.exo_steady_state',ep.periods+2,1);