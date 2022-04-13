function [Y, lhssub, X, residual, fp, lp] = parse_ols_style_equation(ds, ast)
%function [Y, lhssub, X, residual, fp, lp] = parse_ols_style_equation(ds, ast)
% Run OLS on chosen model equations; unlike olseqs, allow for time t
% endogenous variables on LHS
%
% INPUTS
%   ds          [dseries]     data
%   ast         [struct]      AST representing the equation to be parsed
%
% OUTPUTS
%   Y           [dseries]     LHS of the equation (with lhssub subtracted)
%   lhssub      [dseries]     RHS subtracted from LHS
%   X           [dseries]     RHS of the equation
%   residual    [string]      name of residual in equation
%   fp          [date]        first common observed period between Y, lhssub, and X
%   lp          [date]        last common observed period between Y, lhssub, and X
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2019 Dynare Team
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

%% Check inputs
if nargin ~= 2
    error('parse_ols_style_equation takes 2 arguments')
end

if isempty(ds) || ~isdseries(ds)
    error('parse_ols_style_equation: arg 1 must be a dseries');
end

if isempty(ast) || ~isstruct(ast)
    error('parse_ols_style_equation: arg 2 must be a struct');
end

line = ast.line;
if ~strcmp(ast.AST.node_type, 'BinaryOpNode') ...
        || ~strcmp(ast.AST.op, '=')
    parsing_error('expecting equation with equal sign', line);
end

% Check LHS
if ~strcmp(ast.AST.arg1.node_type, 'VariableNode') ...
        && ~strcmp(ast.AST.arg1.node_type, 'UnaryOpNode')
    parsing_error('expecting Variable or UnaryOp on LHS', line);
else
    if ~isOlsVar(ds, ast.AST.arg1) || ~isBaseVarLagEqualToZero(ast.AST.arg1)
        parsing_error('the LHS of the equation must be an Variable or UnaryOp with lag == 0 that exists in the dataset', line);
    end
end

%% Set LHS (Y)
lhssub = dseries();
Y = evalNode(ds, ast.AST.arg1, line, dseries());

%% Set RHS (X)
residual = '';
X = dseries();
terms = decomposeAdditiveTerms([], ast.AST.arg2, 1);
for i = 1:length(terms)
    Xtmp = dseries();
    node_sign = terms{i}{2};
    node_to_parse = terms{i}{1};
    if strcmp(node_to_parse.node_type, 'VariableNode')
        if strcmp(node_to_parse.type, 'parameter')
            % Intercept
            Xtmp = dseries(1, ds.dates, node_to_parse.name)*node_sign;
        elseif strcmp(node_to_parse.type, 'exogenous') && ~any(strcmp(ds.name, node_to_parse.name))
            % Residual if not contained in ds
            if isempty(residual)
                residual = node_to_parse.name;
            else
                parsing_error(['only one residual allowed per equation; encountered ' residual ' & ' node_to_parse.name], line);
            end
        elseif strcmp(node_to_parse.type, 'endogenous') ...
                || (strcmp(node_to_parse.type, 'exogenous') && any(strcmp(ds.name, node_to_parse.name)))
            % Subtract VariableNode from LHS
            % NB: treat exogenous that exist in ds as endogenous
            lhssub = lhssub + evalNode(ds, node_to_parse, line, dseries())*node_sign;
        else
            parsing_error('unexpected variable type found', line, node_to_parse);
        end
    elseif strcmp(node_to_parse.node_type, 'UnaryOpNode')
        % Subtract UnaryOpNode from LHS
        % NB: treat exogenous that exist in ds as endogenous
        lhssub = lhssub + evalNode(ds, node_to_parse, line, dseries());
    elseif strcmp(node_to_parse.node_type, 'BinaryOpNode') && strcmp(node_to_parse.op, '/')
        % Subtract Node from LHS
        % if either arg contains a parameter, it's a parsing error.
        if containsParameter(node_to_parse, line)
            parsing_error('unexpected node found', line, node_to_parse)
        end
        lhssub = lhssub + evalNode(ds, node_to_parse, line, dseries());
    elseif strcmp(node_to_parse.node_type, 'BinaryOpNode') && strcmp(node_to_parse.op, '*')
        % Parse param_expr * endog_expr
        [Xtmp, names] = parseTimesNode(ds, node_to_parse, line);
        Xtmp = Xtmp*node_sign;
        if Xtmp.vobs > 1 || ...
                (Xtmp.vobs == 1 && ~isnan(str2double(Xtmp.name)))
            % Handle constraits
            % Look through Xtmp names for constant
            % if found, subtract from LHS
            for j = 1:length(names)
                if strcmp(names{j}{1}.node_type, 'NumConstNode')
                    pname = num2str(names{j}{1}.value);
                elseif strcmp(names{j}{1}.node_type, 'VariableNode')
                    pname = names{j}{1}.name;
                else
                    parsing_error('unexpected node type', node_to_parse, line);
                end
                psign = names{j}{2};
                if ~isnan(str2double(pname))
                    lhssub = lhssub + psign * str2double(pname) * Xtmp.(pname);
                    Xtmp = Xtmp.remove(pname);
                else
                    % Multiply by psign now so that it can be added together below
                    % Otherwise, it would matter which was encountered first,
                    % a parameter on its own or a linear constraint
                    Xtmp.(pname) = psign * Xtmp.(pname);
                end
            end
        end
    else
        parsing_error('didn''t expect to arrive here', line, node_to_parse);
    end

    names = Xtmp.name;
    for j = length(names):-1:1
        % Handle constraits
        if any(strcmp(X.name, names{j}))
            X.(names{j}) = X.(names{j}) + Xtmp.(names{j});
            Xtmp = Xtmp.remove(names{j});
        end
    end
    X = [X Xtmp];
end
Y = Y - lhssub;

%% Set start and end dates
fp = Y.firstobservedperiod;
lp = Y.lastobservedperiod;
if ~isempty(X)
    % X is empty when AR(1) without parameter is encountered
    fp = max(fp, X.firstobservedperiod);
    lp = min(lp, X.lastobservedperiod);
end
if ~isempty(lhssub)
    fp = max(fp, lhssub.firstobservedperiod);
    lp = min(lp, lhssub.lastobservedperiod);
end

% If it exists, account for tag set in mod file
if isfield(ast, 'tags') ...
        && isfield(ast.tags, 'sample') ...
        && ~isempty(ast.tags.sample)
    colon_idx = strfind(ast.tags.sample, ':');
    fsd = dates(ast.tags.sample(1:colon_idx-1));
    lsd = dates(ast.tags.sample(colon_idx+1:end));
    if fp > fsd
        warning(['The sample over which you want to estimate contains NaNs. '...
            'Adjusting estimation range to begin on: ' fp.char])
    else
        fp = fsd;
    end
    if lp < lsd
        warning(['The sample over which you want to estimate contains NaNs. '...
            'Adjusting estimation range to end on: ' lp.char])
    else
        lp = lsd;
    end
end

Y = Y(fp:lp);
if ~isempty(X)
    X = X(fp:lp);
    names = X.name;
    for i = 1:length(names)
        if all(X.(names{i}).data == 0)
            X = X.remove(names{i});
        end
    end
end
if ~isempty(lhssub)
    lhssub = lhssub(fp:lp);
end
end

%% Helper Functions
function parsing_error(msg, line, node)
if nargin == 3 && ~isempty(node)
    error('\nERROR encountered parsing `%s` in equation on line %d: %s\n', printNode(node), line, msg);
else
    error('\nERROR encountered parsing of equation on line %d: %s\n', line, msg)
end
end

function str = printNode(node)
if strcmp(node.node_type, 'NumConstNode')
    str = num2str(node.value);
elseif strcmp(node.node_type, 'VariableNode')
    if strcmp(node.type, 'endogenous') || strcmp(node.type, 'exogenous')
        str = node.name;
        if node.lag ~= 0
            str = [str '(' num2str(node.lag) ')'];
        end
    elseif strcmp(node.type, 'parameter')
        str = node.name;
    end
elseif strcmp(node.node_type, 'UnaryOpNode')
    str = printNode(node.arg);
    str = [node.op '(' str ')'];
elseif strcmp(node.node_type, 'BinaryOpNode')
    str = ['(' printNode(node.arg1) node.op printNode(node.arg2) ')'];
end
end

function terms = decomposeAdditiveTerms(terms, node, node_sign)
if strcmp(node.node_type, 'NumConstNode') || strcmp(node.node_type, 'VariableNode')
    terms = [terms {{node node_sign}}];
elseif strcmp(node.node_type, 'UnaryOpNode')
    if strcmp(node.op, 'uminus')
        terms = decomposeAdditiveTerms(terms, node.arg, -node_sign);
    else
        terms = [terms {{node node_sign}}];
    end
elseif strcmp(node.node_type, 'BinaryOpNode')
    if strcmp(node.op, '+') || strcmp(node.op, '-')
        terms = decomposeAdditiveTerms(terms, node.arg1, node_sign);
        if strcmp(node.op, '+')
            terms = decomposeAdditiveTerms(terms, node.arg2, node_sign);
        else
            terms = decomposeAdditiveTerms(terms, node.arg2, -node_sign);
        end
    else
        terms = [terms {{node node_sign}}];
    end
else
    terms = [terms {{node node_sign}}];
end
end

function [X, pterms] = parseTimesNode(ds, node, line)
% Separate the parameter expression from the endogenous expression
assert(strcmp(node.node_type, 'BinaryOpNode') && strcmp(node.op, '*'))
if isOlsParamExpr(node.arg1, line)
    pterms = decomposeAdditiveTerms([], node.arg1, 1);
    X = evalNode(ds, node.arg2, line, dseries());
elseif isOlsParamExpr(node.arg2, line)
    pterms = decomposeAdditiveTerms([], node.arg2, 1);
    X = evalNode(ds, node.arg1, line, dseries());
else
    parsing_error('expecting (param expr)*(var expr)', line, node);
end
if strcmp(pterms{1}{1}.node_type, 'NumConstNode')
    X = X.rename(num2str(pterms{1}{1}.value));
elseif strcmp(pterms{1}{1}.node_type, 'VariableNode')
    X = X.rename(pterms{1}{1}.name);
else
    parsing_error('unexpected type', line, node)
end
for ii = 2:length(pterms)
    if strcmp(pterms{ii}{1}.node_type, 'NumConstNode')
        X = [X dseries(X{1}.data, X{1}.firstdate, num2str(pterms{ii}{1}.value))];
    elseif strcmp(pterms{ii}{1}.node_type, 'VariableNode')
        X = [X dseries(X{1}.data, X{1}.firstdate, pterms{ii}{1}.name)];
    else
        parsing_error('unexpected type', line, node)
    end
end
end

function tf = isBaseVarLagEqualToZero(node)
if strcmp(node.node_type, 'VariableNode')
    tf = node.lag == 0;
elseif strcmp(node.node_type, 'UnaryOpNode')
    tf = isBaseVarLagEqualToZero(node.arg);
else
    tf = false;
end
end

function tf = isOlsVar(ds, node)
if strcmp(node.node_type, 'VariableNode') ...
        && (strcmp(node.type, 'endogenous') ...
        || (strcmp(node.type, 'exogenous') && any(strcmp(ds.name, node.name))))
    tf = true;
elseif strcmp(node.node_type, 'UnaryOpNode')
    tf = isOlsVar(ds, node.arg);
else
    tf = false;
end
end

function X = evalNode(ds, node, line, X)
global M_
if strcmp(node.node_type, 'NumConstNode')
    X = dseries(node.value, ds.dates, 'const');
elseif strcmp(node.node_type, 'VariableNode')
    if strcmp(node.type, 'endogenous') ...
            || (strcmp(node.type, 'exogenous') && any(strcmp(ds.name, node.name)))
        if ds.exist(node.name)
            X = ds.(node.name)(node.lag);
        else
            error('Variable %s is not available in the database.', node.name)
        end
    elseif strcmp(node.type, 'parameter')
        X = M_.params(not(cellfun('isempty', strfind(M_.param_names, node.name))));
        if isnan(X) || isinf(X) || ~isreal(X)
            parsing_error(['Value incorrectly set for parameter: ' node.name], line);
        end
    end
elseif strcmp(node.node_type, 'UnaryOpNode')
    Xtmp = evalNode(ds, node.arg, line, X);
    % Only works if dseries supports . notation for unary op (true for log/diff)
    % Otherwise, use: X = eval([node.op '(Xtmp)']);
    try
        if strcmp(node.op, 'uminus')
            X = -Xtmp;
        else
            X = Xtmp.(node.op);
        end
        if any(isinf(X)) || ~isreal(X)
            parsing_error(['Error applying ' node.op], line, node);
        end
    catch
        parsing_error(['Error applying ' node.op], line, node);
    end
elseif strcmp(node.node_type, 'BinaryOpNode')
    Xtmp1 = evalNode(ds, node.arg1, line, X);
    Xtmp2 = evalNode(ds, node.arg2, line, X);
    switch node.op
        case '*'
            Xtmp = Xtmp1 * Xtmp2;
        case '/'
            Xtmp = Xtmp1 / Xtmp2;
        case '+'
            Xtmp = Xtmp1 + Xtmp2;
        case '-'
            Xtmp = Xtmp1 - Xtmp2;
        otherwise
            parsing_error(['got unexpected binary op ' node.op], line, node);
    end
    if any(isinf(Xtmp)) || ~isreal(Xtmp)
        parsing_error(['Error applying ' node.op], line, node);
    end
    X = X + Xtmp;
else
    parsing_error(['got unexpected node type ' node.node_type], line, node);
end
end

function tf = isOlsParamExpr(node, line)
if strcmp(node.node_type, 'NumConstNode')
    tf = true;
elseif strcmp(node.node_type, 'VariableNode')
    if strcmp(node.type, 'parameter')
        tf = true;
    else
        tf = false;
    end
elseif strcmp(node.node_type, 'UnaryOpNode')
    tf = false;
elseif strcmp(node.node_type, 'BinaryOpNode')
    tf = isOlsParamExpr(node.arg1, line) && isOlsParamExpr(node.arg2, line);
    if tf && ~strcmp(node.op, '-')
        parsing_error(['got unexpected op ' node.op], line, node);
    end
else
    parsing_error(['got unexpected type ' node.node_type], line, node);
end
end

function tf = containsParameter(node, line)
if strcmp(node.node_type, 'NumConstNode')
    tf = false;
elseif strcmp(node.node_type, 'VariableNode')
    if strcmp(node.type, 'parameter')
        tf = true;
    else
        tf = false;
    end
elseif strcmp(node.node_type, 'UnaryOpNode')
    tf = containsParameter(node.arg, line);
elseif strcmp(node.node_type, 'BinaryOpNode')
    tf = containsParameter(node.arg1, line) || containsParameter(node.arg2, line);
else
    parsing_error(['got unexpected type ' node.node_type], line, node);
end
end
