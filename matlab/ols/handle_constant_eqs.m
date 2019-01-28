function [ast, ds] = handle_constant_eqs(ast, ds)
%function [ast, ds] = handle_constant_eqs(ast, ds)
%
% Code to handle equations of type `X = 0;` in AST. Returns equation(s)
% removed from AST and ds.X == 0.
%
% INPUTS
%   ast       [cell array]  JSON representation of abstract syntax tree
%   ds        [dseries]     data to be updated
%
% OUTPUTS
%   ast       [cell array]  updated JSON representation of abstract syntax tree
%   ds        [dseries]     data to be updated
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2019 Dynare Team
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

if nargin ~= 2
    error('Incorrect number of arguments to function')
end

if isempty(ast) || ~iscell(ast)
    error('ast must be a cell')
end

if isempty(ds) || ~isdseries(ds)
    error('ds must be a nonempty dseries')
end

for i = length(ast):-1:1
    assert(strcmp(ast{i}.AST.node_type, 'BinaryOpNode') && strcmp(ast{i}.AST.op, '='), 'Expecting equation');
    if strcmp(ast{i}.AST.arg2.node_type, 'NumConstNode')
        if ~strcmp(ast{i}.AST.arg1.node_type, 'VariableNode')
            error('At the moment only handling Variable Nodes on LHS')
        end
        ds.(ast{i}.AST.arg1.name) = ...
            dseries(ast{i}.AST.arg2.value*ones(ds.nobs, 1), ds.firstdate, ast{i}.AST.arg1.name);
        ast = ast([1:i-1, i+1:end]);
    end
end
end
