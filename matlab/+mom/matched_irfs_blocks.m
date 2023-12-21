function [data_irfs, weight_mat, irf_index, max_irf_horizon] = matched_irfs_blocks(matched_irfs, matched_irfs_weight, varobs_id, obs_nbr, exo_nbr, endo_names, exo_names)
% [data_irfs, weight_mat, irf_index, max_irf_horizon] = matched_irfs_blocks(matched_irfs, matched_irfs_weight, varobs_id, obs_nbr, exo_nbr, endo_names, exo_names)
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
% endo_names:          [cell array] names of endogenous variables
% exo_names:           [cell array] names of exogenous variables
% -------------------------------------------------------------------------
% OUTPUT
% data_irfs:           [matrix]     IRFs for VAROBS as declared in matched_irfs block
% weight_mat:          [matrix]     weighting matrix for IRFs as declared in matched_irfs_weight block
% irf_index:           [vector]     index for selecting specific IRFs from full IRF matrix of observables
% max_irf_horizon:     [scalar]     maximum IRF horizon as declared in matched_irfs block
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

% note matched_irfs block:
% - each row in the cell contains a unique combination of var and varexo,
%   however the third column in each row is a nested cell with information
%   on periods, values and weights
% - periods, values and weights can span several rows with different lengths of entries
% - in some cases we need to duplicate values and/or weights
% - at the end we want to have everything vectorized and the same length

% get maximum IRF horizons
max_irf_horizon = [];
for jj = 1:size(matched_irfs,1)
    max_irf_horizon = [max_irf_horizon; cell2mat(cellfun(@(c) c(:), matched_irfs{jj,3}(:,1), 'UniformOutput', false))];
end
max_irf_horizon = max(max_irf_horizon);

% create full matrix where 1st dimension are IRF periods, 2nd dimension are variables as declared in VAROBS, 3rd dimension are shocks
% idea: overwrite NaN values if they are declared in matched_irfs block; at the end the remaining NaN values will be removed
data_irfs = NaN(max_irf_horizon,obs_nbr,exo_nbr);
% create full empirical weighting matrix, identity matrix by default, i.e. all IRFs are equally important
% idea: first specify full matrix and then reduce it using only entries that are declared in matched_irfs block
weight_mat = speye(max_irf_horizon*obs_nbr*exo_nbr);

for jj = 1:size(matched_irfs,1)
    id_var       = find(ismember(endo_names,matched_irfs{jj,1}));
    id_varobs    = find(varobs_id==id_var,1);
    id_shock = find(ismember(exo_names,matched_irfs{jj,2}));
    if isempty(id_varobs)
        skipline;
        error('method_of_moments: You specified an IRF matching involving variable %s, but it is not declared as a varobs!',endo_names{id_var})
    end    
    IRF_PERIODS = []; IRF_VALUES = []; IRF_WEIGHTS = [];
    for kk = 1:size(matched_irfs{jj,3},1)
        irf_periods = matched_irfs{jj,3}{kk,1};
        if length(unique(irf_periods)) < length(irf_periods) % row-specific check for unique periods
            error('method_of_moments: You specified an IRF matching involving variable %s and shock %s, but there were duplicate ''periods'' in the specification!',endo_names{id_var},exo_names{id_shock});
        end
        irf_values = matched_irfs{jj,3}{kk,2};
        if length(irf_values)==1
            irf_values = repmat(irf_values,length(irf_periods),1);
        end
        if length(irf_periods) ~= length(irf_values) % row-specific check for enough values
            error('method_of_moments: You specified an IRF matching involving variable %s and shock %s, but the length of ''periods'' does not match the length of ''values''!',endo_names{id_var},exo_names{id_shock});
        end
        irf_weights = matched_irfs{jj,3}{kk,3};
        if length(irf_weights)==1
            irf_weights = repmat(irf_weights,length(irf_periods),1);
        end
        if length(irf_periods) ~= length(irf_weights) % row-specific check for enough weights
            error('method_of_moments: You specified an IRF matching involving variable %s and shock %s, but the length of ''periods'' does not match the length of ''weights''!',endo_names{id_var},exo_names{id_shock});
        end
        IRF_PERIODS = [IRF_PERIODS; irf_periods(:)];
        IRF_VALUES = [IRF_VALUES; irf_values(:)];
        IRF_WEIGHTS = [IRF_WEIGHTS; irf_weights(:)];
    end
    if length(unique(irf_periods)) < length(irf_periods) % overall check for unique periods
        error('method_of_moments: You specified an IRF matching involving variable %s and shock %s, but there were duplicate ''periods'' in the specification!',endo_names{id_var},exo_names{id_shock});
    end
    for hh = 1:length(IRF_PERIODS)
        data_irfs(IRF_PERIODS(hh),id_varobs,id_shock) = IRF_VALUES(hh);
        if IRF_WEIGHTS(hh) ~= 1
            idweight_mat = sub2ind(size(data_irfs),IRF_PERIODS(hh),id_varobs,id_shock);
            weight_mat(idweight_mat,idweight_mat) = IRF_WEIGHTS(hh);
        end
    end
end

% fine-tune weighting matrix using matched_irfs_weights
for jj = 1:size(matched_irfs_weight,1)
    id_var1 = find(ismember(endo_names,matched_irfs_weight{jj,1}));
    id_var2 = find(ismember(endo_names,matched_irfs_weight{jj,4}));
    id_varobs1 = find(varobs_id==id_var1,1);
    id_varobs2 = find(varobs_id==id_var2,1);
    if isempty(id_varobs1)
        skipline;
        error('method_of_moments: You specified a weight for an IRF matching involving variable %s, but it is not a varobs!',endo_names{id_var1})
    end
    if isempty(id_varobs2)
        skipline;
        error('method_of_moments: You specified a weight for an IRF matching involving variable %s, but it is not a varobs!',endo_names{id_var2})
    end
    id_shock1 = find(ismember(exo_names,matched_irfs_weight{jj,3}));
    id_shock2 = find(ismember(exo_names,matched_irfs_weight{jj,6}));
    irf_periods1 = matched_irfs_weight{jj,2};
    irf_periods2 = matched_irfs_weight{jj,5};
    if length(irf_periods1) ~= length(irf_periods2)
        error('method_of_moments: You specified a ''matched_irfs_weights'' entry for an IRF matching involving %s/%s and %s/%s,\n                   but the horizons do not have the same length!',endo_names{id_var1},exo_names{id_shock1},endo_names{id_var2},exo_names{id_shock2});
    end
    if max([irf_periods1(:);irf_periods2(:)]) > max_irf_horizon
        error('method_of_moments: You specified a ''matched_irfs_weights'' entry for an IRF matching involving %s/%s and %s/%s,\n                   but the horizon is larger than the maximum one declared in the ''matched_irfs'' block!',endo_names{id_var1},exo_names{id_shock1},endo_names{id_var2},exo_names{id_shock2});
    end    
    weight_mat_values = matched_irfs_weight{jj,7};
    if length(weight_mat_values)==1 && length(irf_periods1)>1
        weight_mat_values = repmat(weight_mat_values,length(irf_periods1),1);
    end
    if length(weight_mat_values) ~= length(irf_periods1)
        error('method_of_moments: You specified a ''matched_irfs_weights'' entry for an IRF matching involving %s/%s and %s/%s,\n                   but the horizons do not match the length of ''weights''!',endo_names{id_var1},exo_names{id_shock1},endo_names{id_var2},exo_names{id_shock2});
    end
    for hh = 1:length(irf_periods1)
        idweight_mat1 = sub2ind(size(data_irfs),irf_periods1(hh),id_varobs1,id_shock1);
        idweight_mat2 = sub2ind(size(data_irfs),irf_periods2(hh),id_varobs2,id_shock2);
        weight_mat(idweight_mat1,idweight_mat2) = weight_mat_values(hh);
        weight_mat(idweight_mat2,idweight_mat1) = weight_mat_values(hh); % symmetry
    end
end

% remove non-specified IRFs
irf_index = find(~isnan(data_irfs));
data_irfs = data_irfs(irf_index);
weight_mat = weight_mat(irf_index,irf_index);