function expression = remove_aux_variables_from_expression(expression, M_)
% expression = remove_aux_variables_from_expression(expression, M_)

% Copyright © 2022-2023 Dynare Team
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

% Get list of endogenous variables in expression
list_of_words = regexp(expression, '\<\w*\>', 'match');
list_of_words = setdiff(list_of_words, M_.param_names);
list_of_words = setdiff(list_of_words, M_.exo_names);
isnotanumber = isnan(str2double(list_of_words));
list_of_words = list_of_words(isnotanumber);
list_of_words = setdiff(list_of_words, {'diff','log','exp'});

for i=1:length(list_of_words)
    id = find(strcmp(list_of_words{i}, M_.endo_names));
    if isempty(id) || id<=M_.orig_endo_nbr
        continue
    end
    auxinfo = M_.aux_vars(get_aux_variable_id(id));
    expression = regexprep(expression, sprintf('\\<%s\\>', list_of_words{i}), auxinfo.orig_expr);
end