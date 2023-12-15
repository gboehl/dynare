function indx = get_new_or_existing_ei_index(substructure_name, name1, name2)
% function indx = get_new_or_existing_ei_index(substructure_name, name1, name2)
%
% Returns the new estimation_info.substructure_name index
% for the name1 & name2 pair
%
% INPUTS
%    substructure_name  [string]  the name of the substructure in which to look
%    name1              [string]  the variable for which the subsample was declared
%    name2              [string]  used only in case of corr(name1,name2).subsamples()
%
% OUTPUTS
%    indx               [integer] new index in the
%                                 estimation_info.substructure_name
%                                 associated with the name pair
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2012-2017 Dynare Team
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

global estimation_info

if isempty(estimation_info.(substructure_name))
    indx = 1;
    return
end

if isempty(name2) % parameter or std() statement
    indx = find(strcmp(name1, estimation_info.(substructure_name)));
else % corr statement
    indx = find(strcmp(sprintf('%s:%s', name1, name2), estimation_info.(substructure_name)));
    if isempty(indx)
        indx = find(strcmp(sprintf('%s:%s', name2, name1), estimation_info.(substructure_name)));
    end
end

if isempty(indx)
    indx = size(estimation_info.(substructure_name), 2) + 1;
end

if size(indx, 2) > 1
    error('Error: %s %s found more than once in estimation_info.subsamples_index', name1, name2);
end
