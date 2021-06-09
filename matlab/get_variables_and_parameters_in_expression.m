function objects = get_variables_and_parameters_in_expression(expr)

% Returns the variables and parameters appearing in an expression.
%
% INPUTS
% - expr       [char]             1×m char array, dynare model expression (typically RHS or LHS of an equation).
%
% OUTPUTS
% - objects    [cell]             cell of row char arrays, names of the variables and parameters in expr.

% Copyright (C) 2020 Dynare Team
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

objects = strsplit(expr, {'+','-','*','/','^', ...
                    'log(', 'log10(', 'ln(', 'exp(', ...
                    'sqrt(', 'abs(', 'sign(', ...
                    'sin(', 'cos(', 'tan(', 'asin(', 'acos(', 'atan(', ...
                    'min(', 'max(', ...
                    'normcdf(', 'normpdf(', 'erf(', ...
                    'diff(', 'adl(', '(', ')', '\n', '\t', ' '});

% Filter out the numbers, punctuation.
objects(cellfun(@(x) all(isstrprop(x, 'digit')+isstrprop(x, 'punct')), objects)) = [];

% Filter out empty elements.
objects(cellfun(@(x) all(isempty(x)), objects)) = [];