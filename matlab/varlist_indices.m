function [i_var, nvar, index_unique_present] = varlist_indices(sublist, list, nocheck_dummy)

% returns the indices of a list of variables
%
% INPUT
%   sublist     [cell of char arrays] sublist of variables
%   list        [cell of char arrays] list of variables
%
% OUTPUT
%   i_var                             variable indices in list
%   nvar                              number of unique variables contained in list
%   index_uniques                     indices of unique and present elements in sublist
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2010-2021 Dynare Team
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

if nargin<3
    nocheck_dummy=0;
end
if isempty(sublist)
    check = [];
    i_var = [];
else
    [check, i_var] = ismember(sublist, list);
end

indices_not_present=[];
if ~all(check)
    k_not_present = find(~check);
    if ~nocheck_dummy
        str = 'The following symbols are not endogenous variables:';
        for ii = 1:length(k_not_present)
            str = sprintf('%s %s', str, sublist{k_not_present(ii)});
        end
        error(str)
    else
        indices_not_present=find(i_var==0);
    end
end

nvar_present = length(i_var(check));
[~, index_unique] = unique(i_var, 'first');
index_unique_present = index_unique(~ismember(index_unique,indices_not_present));
index_unique_present = sort(index_unique_present);
i_var_unique_present = i_var(index_unique_present);

if length(i_var_unique_present)~=nvar_present
    k = find(~ismember((1:length(i_var))',index_unique_present) & i_var~=0);
    str = 'The following symbols are specified twice in the variable list and are considered only once:';
    for ii = 1:length(k)
        str = sprintf('%s %s', str, sublist{k(ii)});
    end
    warning('%s\n', str)
end

i_var = i_var_unique_present;
nvar = length(i_var);
