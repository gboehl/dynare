function matched_moments = matched_moments_block(matched_moments, mom_method)
% matched_moments = matched_moments_block(matched_moments, mom_method)
% -------------------------------------------------------------------------
% Checks and transforms matched_moments block for further use in the estimation
% -------------------------------------------------------------------------
% INPUTS
% matched_moments: [cell array] original matched_moments block
% mom_method:      [string]     method of moments method (GMM or SMM)
% -------------------------------------------------------------------------
% OUTPUT
% matched_moments: [cell array] transformed matched_moments block
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


matched_moments_orig = matched_moments;
% higher-order product moments not supported yet for GMM
if strcmp(mom_method, 'GMM') && any(cellfun(@sum,matched_moments(:,3))> 2)
    error('method_of_moments: GMM does not yet support product moments higher than 2. Change your ''matched_moments'' block!');
end
% check for duplicate moment conditions
for jm = 1:size(matched_moments,1)
    % expand powers to vector of ones
    if any(matched_moments{jm,3}>1)
        tmp1=[]; tmp2=[]; tmp3=[];
        for jjm=1:length(matched_moments{jm,3})
            tmp1 = [tmp1 repmat(matched_moments{jm,1}(jjm),[1 matched_moments{jm,3}(jjm)]) ];
            tmp2 = [tmp2 repmat(matched_moments{jm,2}(jjm),[1 matched_moments{jm,3}(jjm)]) ];
            tmp3 = [tmp3 repmat(1,[1 matched_moments{jm,3}(jjm)]) ];
        end
        matched_moments{jm,1} = tmp1;
        matched_moments{jm,2} = tmp2;
        matched_moments{jm,3} = tmp3;
    end
    % shift time structure to focus only on lags
    matched_moments{jm,2} = matched_moments{jm,2} - max(matched_moments{jm,2});
    % sort such that t=0 variable comes first
    [matched_moments{jm,2},idx_sort] = sort(matched_moments{jm,2},'descend');
    matched_moments{jm,1} = matched_moments{jm,1}(idx_sort);
    matched_moments{jm,3} = matched_moments{jm,3}(idx_sort);
end
% find duplicate rows in cell array by making groups according to powers as we can then use cell2mat for the unique function
powers = cellfun(@sum,matched_moments(:,3))';
unique_mom_idx = [];
for jpow = unique(powers)
    idx1 = find(powers==jpow);
    [~,idx2] = unique(cell2mat(matched_moments(idx1,:)),'rows');
    unique_mom_idx = [unique_mom_idx idx1(idx2)];
end
% remove duplicate elements
duplicate_moms = setdiff(1:size(matched_moments_orig,1),unique_mom_idx);
if ~isempty(duplicate_moms)
    fprintf('Duplicate declared moments found and removed in ''matched_moments'' block in rows:\n %s.\n',num2str(duplicate_moms));
    fprintf('Dynare will continue with remaining moment conditions\n');
end
if strcmp(mom_method, 'SMM')
    % for SMM: keep the original structure, but get rid of duplicate moments
    matched_moments = matched_moments_orig(sort(unique_mom_idx),:);
elseif strcmp(mom_method, 'GMM')
    % for GMM we use the transformed matched_moments structure
    matched_moments = matched_moments(sort(unique_mom_idx),:);
end