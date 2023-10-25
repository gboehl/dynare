function [expression, lhs] = rewrite_equation_with_tables(expression, lhs, islaggedvariables, pnames, enames, xnames, pid, eid, xid, data)

% Takes and expression of parameters and variables (with max. lag equal to one)
% and replace parameters by elements of a vector and variables by columns of a
% data matrix.

% Copyright © 2021 Dynare Team
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

% List of objects to be replaced
objNames = [pnames; enames; xnames];
objIndex = [pid; eid; xid];
objTypes = [ones(length(pid), 1); 2*ones(length(eid)+length(xid), 1)];

[~,I] = sort(cellfun(@length, objNames), 'descend');
objNames = objNames(I);
objIndex = objIndex(I);
objTypes = objTypes(I);

% Substitute parameters and variables.  Note that in the transformed model (we consider rhs
% instead of RHS) the maximum lag is 1.
for i=1:length(objNames)
    switch objTypes(i)
      case 1
        expression = strrep(expression, objNames{i}, sprintf('M_.params(%u)', objIndex(i)));
      case 2
        k = find(strcmp(objNames{i}, data.name));
        if isempty(k)
            error('Variable %s is missing in the database.', objNames{i})
        end
        j = regexp(expression, ['\<', objNames{i}, '\>']);
        if islaggedvariables
            jlag = regexp(expression, ['\<', objNames{i}, '\(-1\)']);
            if ~isempty(jlag)
                expression = regexprep(expression, ['\<' objNames{i} '\(-1\)'], sprintf('data(1:end-1,%u)', k));
            end
            if ~isempty(setdiff(j, jlag))
                expression = regexprep(expression, ['\<' objNames{i} '\>'], sprintf('data(2:end,%u)', k));
            end
        else
            expression = regexprep(expression, ['\<' objNames{i} '\>'], sprintf('data(:,%u)', k));
        end
        if contains(lhs, objNames{i})
            if islaggedvariables
                lhs = strrep(lhs, objNames{i}, sprintf('data(2:end,%u)', k));
            else
                lhs = strrep(lhs, objNames{i}, sprintf('data(:,%u)', k));
            end
        end
    end
end

% Allow elementwise operations
expression = strrep(expression, '^', '.^');
expression = strrep(expression, '/', './');
expression = strrep(expression, '*', '.*');