function [data_irfs, weight_mat, irf_index, max_irf_horizon] = matched_irfs_blocks(matched_irfs, matched_irfs_weight, varobs_id, obs_nbr, exo_nbr, endo_names)
% [data_irfs, weight_mat, irf_index, max_irf_horizon] = matched_irfs_blocks(matched_irfs, matched_irfs_weight, varobs_id, obs_nbr, exo_nbr, endo_names)
% -------------------------------------------------------------------------
% Checks and transforms matched_irfs and matched_irfs_weight blocks
% for further use in the estimation.
% -------------------------------------------------------------------------
% INPUTS
% matched_irfs:        [cell array] original matched_irfs block
% matched_irfs_weight: [cell array] original matched_irfs_weight block
% varobs_id:           [vector]     index for observable variables in endo_names
% obs_nbr:             [scalar]     number of observable variables
% exo_nbr:             [scalar]     number of exogenous variables
% endo_names:          [cell array] list of endogenous variables
% -------------------------------------------------------------------------
% OUTPUT
% data_irfs:           [matrix]     irfs for VAROBS as declared in matched_irfs block
% weight_mat:          [matrix]     weighting matrix for irfs as declared in matched_irfs_weight block
% irf_index:           [vector]     index for selecting specific irfs from full irf matrix of observables
% max_irf_horizon:     [scalar]     maximum irf horizon as declared in matched_irfs block
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


max_irf_horizon = max(cellfun(@(x) x(end), matched_irfs(:,1))); % get maximum irf horizon
% create full matrix where 1st dimension are irf periods, 2nd dimension are variables as declared in VAROBS, 3rd dimension are shocks.
data_irfs = nan(max_irf_horizon,obs_nbr,exo_nbr);
% overwrite nan values if they are declared in matched_irfs block; remaining nan values will be later ignored in the matching
for jj = 1:size(matched_irfs,1)
    id_var       = matched_irfs{jj,1}(1);
    id_varobs    = find(varobs_id==id_var,1);
    id_shock     = matched_irfs{jj,1}(2);
    id_irf_period = matched_irfs{jj,1}(3);
    irf_value    = matched_irfs{jj,2};
    if isempty(id_varobs)
        skipline;
        error('method_of_moments: You specified an irf matching involving variable %s, but it is not declared as a varobs!',endo_names{id_var})
    end
    data_irfs(id_irf_period,id_varobs,id_shock) = irf_value;
end
% create (full) empirical weighting matrix
weight_mat = eye(max_irf_horizon*obs_nbr*exo_nbr); % identity matrix by default: all irfs are equally important
for jj = 1:size(matched_irfs_weight,1)
    id_var1 = matched_irfs_weight{jj,1}(1);  id_varobs1 = find(varobs_id==id_var1,1);  id_shock1 = matched_irfs_weight{jj,1}(2);  id_irf_period1 = matched_irfs_weight{jj,1}(3);
    id_var2 = matched_irfs_weight{jj,2}(1);  id_varobs2 = find(varobs_id==id_var2,1);  id_shock2 = matched_irfs_weight{jj,2}(2);  id_irf_period2 = matched_irfs_weight{jj,2}(3);
    weight_mat_value = matched_irfs_weight{jj,3};
    if isempty(id_varobs1)
        skipline;
        error('method_of_moments: You specified a weight for an irf matching involving variable %s, but it is not a varobs!',endo_names{id_var1})
    end
    if isempty(id_varobs2)
        skipline;
        error('method_of_moments: You specified a weight for an irf matching involving variable %s, but it is not a varobs!',endo_names{id_var2})
    end
    idweight_mat1 = sub2ind(size(data_irfs),id_irf_period1,id_varobs1,id_shock1);
    idweight_mat2 = sub2ind(size(data_irfs),id_irf_period2,id_varobs2,id_shock2);
    weight_mat(idweight_mat1,idweight_mat2) = weight_mat_value;
    weight_mat(idweight_mat2,idweight_mat1) = weight_mat_value; % symmetry
end
% focus only on specified irfs
irf_index = find(~isnan(data_irfs));
data_irfs = data_irfs(irf_index);
weight_mat = weight_mat(irf_index,irf_index);