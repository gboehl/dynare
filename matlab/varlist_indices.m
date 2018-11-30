function [i_var, nvar, index_uniques] = varlist_indices(sublist, list)

% returns the indices of a list of endogenous variables
%
% INPUT
%   sublist     [cell of char arrays] sublist of variables
%   list        [cell of char arrays] list of variables
%
% OUTPUT
%   i_var                             variable indices in M_.endo_names
%   nvar                              number of variables in varlist
%   index_uniques                     indices of unique elements in varlist
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2010-2018 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if isempty(sublist)
    check = [];
    i_var = [];
else
    [check, i_var] = ismember(sublist, list);
end

if ~all(check)
    k = find(~check);
    str = 'The following symbols are not endogenous variables:';
    for ii = 1:length(k)
        str = sprintf('%s %s', str, sublist{k(ii)});
    end
    error(str)
end

nvar = length(i_var);
[i_var_unique, index_uniques, ~] = unique(i_var, 'first');
index_uniques = sort(index_uniques);
i_var_unique = i_var(index_uniques);

if length(i_var_unique)~=nvar
    k = find(~ismember(1:nvar,index_uniques));
    str = 'The following symbols are specified twice in the variable list and are considered only once:';
    for ii = 1:length(k)
        str = sprintf('%s %s', str, sublist{k(ii)});
    end
    warning('%s\n', str)
    i_var = i_var_unique;
    nvar = length(i_var);
end