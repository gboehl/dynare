function [ast, jsonmodel, ds] = handle_constant_eqs(ast, jsonmodel, ds)
%function [ast, jsonmodel, ds] = handle_constant_eqs(ast, jsonmodel, ds)
%
% Code to handle equations of type `X = 0;` in AST. Returns equation(s)
% removed from AST and ds.X == 0.
%
% INPUTS
%   ast                  [cell array]  JSON representation of abstract syntax tree
%   jsonmodel            [cell array]  JSON representation of model block
%   ds                   [dseries]     data to be updated
%
% OUTPUTS
%   ast                  [cell array]  updated JSON representation of abstract syntax tree
%   jsonmodel            [cell array]  JSON representation of model block
%   ds                   [dseries]     data to be updated
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

for i = length(ast):-1:1
    assert(strcmp(ast{i}.AST.node_type, 'BinaryOpNode') && strcmp(ast{i}.AST.op, '='), 'Expecting equation');
    if strcmp(ast{i}.AST.arg2.node_type, 'NumConstNode')
        if ~strcmp(ast{i}.AST.arg1.node_type, 'VariableNode')
            error('At the moment only handling Variable Nodes on LHS')
        end
        ds.(ast{i}.AST.arg1.name) = ...
            dseries(ast{i}.AST.arg2.value*ones(ds.nobs, 1), ds.firstdate, ast{i}.AST.arg1.name);
        ast = ast([1:i-1, i+1:end]);
        jsonmodel = jsonmodel([1:i-1, i+1:end]);
    end
end
% for i = 1:length(ast)
%     ast{i}.AST.arg2 = replace_in_equtaion_recursive(ast{i}.AST.arg2, vars_and_vals);
% end
end

% function node = replace_in_equtaion_recursive(node, vars_and_vals)
% if strcmp(node.node_type, 'NumConstNode')
% elseif strcmp(node.node_type, 'VariableNode')
%     if strcmp(node.type, 'endogenous')
%         if any(strcmp(vars_and_vals(:,1), node.name))
%             node = struct('node_type', 'NumConstNode', ...
%                 'value', vars_and_vals{strcmp(vars_and_vals(:,1), node.name), 2});
%         end
%     end
% elseif strcmp(node.node_type, 'UnaryOpNode')
%     node.arg = replace_in_equtaion_recursive(node.arg, vars_and_vals);
%     if strcmp(node.arg.node_type, 'NumConstNode')
%         switch node.op
%             case 'log'
%                 node = struct('node_type', 'NumConstNode', ...
%                     'value', log(node.arg.value));
%             case 'exp'
%                 node = struct('node_type', 'NumConstNode', ...
%                     'value', exp(node.arg.value));
%             case 'uminus'
%                 node = struct('node_type', 'NumConstNode', ...
%                     'value', -node.arg.value);
%         end
%     end
% elseif strcmp(node.node_type, 'BinaryOpNode')
%     node.arg1 = replace_in_equtaion_recursive(node.arg1, vars_and_vals);
%     node.arg2 = replace_in_equtaion_recursive(node.arg2, vars_and_vals);
%     if strcmp(node.arg1.node_type, 'NumConstNode')
%         if strcmp(node.arg2.node_type, 'NumConstNode')
%             node = struct('node_type', 'NumConstNode', ...
%                 'value', eval(['node.arg1.value' node.op 'node.arg1.value']));
%         elseif node.arg1.value == 0
%             if strcmp(node.op, 'minus')
%                 node = struct('node_type', 'UnaryOpNode', ...
%                     'op', 'uminus', ...
%                     'arg', node.arg2);
%             elseif strcmp(node.op, 'plus')
%                 node = node.arg2;
%             else
%                 error(['Node type not supported ' node.op]);
%             end
%         end
%     end
%     if strcmp(node.arg2.node_type, 'NumConstNode') && node.arg2.value == 0
%         node = node.arg1;
%     end
% else
%     error(['Node type not supported: ' node.node_type])
% end
% end
