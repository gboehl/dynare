function [data_irfs, weightMat, irfIndex, maxIrfHorizon] = matched_irfs_blocks(matched_irfs, matched_irfs_weight, varobs_id, obs_nbr, exo_nbr, endo_names)
% [data_irfs, weightMat, irfIndex, maxIrfHorizon] = matched_irfs_blocks(matched_irfs, matched_irfs_weight, varobs_id, obs_nbr, exo_nbr, endo_names)
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
% weightMat:           [matrix]     weighting matrix for irfs as declared in matched_irfs_weight block
% irfIndex:            [vector]     index for selecting specific irfs from full irf matrix of observables
% maxIrfHorizon:       [scalar]     maximum irf horizon as declared in matched_irfs block
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


maxIrfHorizon = max(cellfun(@(x) x(end), matched_irfs(:,1))); % get maximum irf horizon
% create full matrix where 1st dimension are irf periods, 2nd dimension are variables as declared in VAROBS, 3rd dimension are shocks.
data_irfs = nan(maxIrfHorizon,obs_nbr,exo_nbr);
% overwrite nan values if they are declared in matched_irfs block; remaining nan values will be later ignored in the matching
for jj = 1:size(matched_irfs,1)
    idVar       = matched_irfs{jj,1}(1);
    idVarobs    = find(varobs_id==idVar,1);
    idShock     = matched_irfs{jj,1}(2);
    idIrfPeriod = matched_irfs{jj,1}(3);
    irfValue    = matched_irfs{jj,2};
    if isempty(idVarobs)
        skipline;
        error('method_of_moments: You specified an irf matching involving variable %s, but it is not declared as a varobs!',endo_names{idVar})
    end
    data_irfs(idIrfPeriod,idVarobs,idShock) = irfValue;
end
% create (full) empirical weighting matrix
weightMat = eye(maxIrfHorizon*obs_nbr*exo_nbr); % identity matrix by default: all irfs are equally important
for jj = 1:size(matched_irfs_weight,1)
    idVar1 = matched_irfs_weight{jj,1}(1);  idVarobs1 = find(varobs_id==idVar1,1);  idShock1 = matched_irfs_weight{jj,1}(2);  idIrfPeriod1 = matched_irfs_weight{jj,1}(3);
    idVar2 = matched_irfs_weight{jj,2}(1);  idVarobs2 = find(varobs_id==idVar2,1);  idShock2 = matched_irfs_weight{jj,2}(2);  idIrfPeriod2 = matched_irfs_weight{jj,2}(3);
    weightMatValue = matched_irfs_weight{jj,3};
    if isempty(idVarobs1)
        skipline;
        error('method_of_moments: You specified a weight for an irf matching involving variable %s, but it is not a varobs!',endo_names{idVar1})
    end
    if isempty(idVarobs2)s
        skipline;
        error('method_of_moments: You specified a weight for an irf matching involving variable %s, but it is not a varobs!',endo_names{idVar2})
    end
    idweightMat1 = sub2ind(size(data_irfs),idIrfPeriod1,idVarobs1,idShock1);
    idweightMat2 = sub2ind(size(data_irfs),idIrfPeriod2,idVarobs2,idShock2);
    weightMat(idweightMat1,idweightMat2) = weightMatValue;
    weightMat(idweightMat2,idweightMat1) = weightMatValue; % symmetry
end
% focus only on specified irfs
irfIndex = find(~isnan(data_irfs));
data_irfs = data_irfs(irfIndex);
weightMat = weightMat(irfIndex,irfIndex);