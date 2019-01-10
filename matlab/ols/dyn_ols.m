function varargout = dyn_ols(ds, fitted_names_dict, eqtags)
% function varargout = dyn_ols(ds, fitted_names_dict, eqtags)
% Run OLS on chosen model equations; unlike olseqs, allow for time t
% endogenous variables on LHS
%
% INPUTS
%   ds                [dseries]    data
%   fitted_names_dict [cell]       Nx2 or Nx3 cell array to be used in naming fitted
%                                  values; first column is the equation tag,
%                                  second column is the name of the
%                                  associated fitted value, third column
%                                  (if it exists) is the function name of
%                                  the transformation to perform on the
%                                  fitted value.
%   eqtags            [cellstr]    names of equation tags to estimate. If empty,
%                                  estimate all equations
%
% OUTPUTS
%   varargout{1}      [dseries]    data updated with fitted values (if not
%                                  called from olsgibbs)
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2017-2018 Dynare Team
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

global M_ oo_ options_

assert(nargin >= 1 && nargin <= 3, 'dyn_ols: takes between 1 and 3 arguments');
assert(isdseries(ds), 'dyn_ols: the first argument must be a dseries');

jsonfile = [M_.fname filesep() 'model' filesep() 'json' filesep() 'modfile-original.json'];
if exist(jsonfile, 'file') ~= 2
    error('Could not find %s! Please use the json=compute option (See the Dynare invocation section in the reference manual).', jsonfile);
end

%% Get Equation(s)
jsonmodel = loadjson(jsonfile);
ast       = jsonmodel.abstract_syntax_tree;
jsonmodel = jsonmodel.model;

if nargin == 1
    fitted_names_dict = {};
elseif nargin == 2
    assert(isempty(fitted_names_dict) || ...
        (iscell(fitted_names_dict) && ...
        (size(fitted_names_dict, 2) == 2 || size(fitted_names_dict, 2) == 3)), ...
        'dyn_ols: the second argument must be an Nx2 or Nx3 cell array');
elseif nargin == 3
    ast = getEquationsByTags(ast, 'name', eqtags);
    jsonmodel = getEquationsByTags(jsonmodel, 'name', eqtags);
end

%% Estimation

for i = 1:length(ast)
    if ~strcmp(ast{i}.AST.node_type, 'BinaryOpNode') ...
            || ~strcmp(ast{i}.AST.op, '=')
        ols_error('expecting equation with equal sign', line);
    end

    % Check LHS
    if ~strcmp(ast{i}.AST.arg1.node_type, 'VariableNode') ...
            && ~strcmp(ast{i}.AST.arg1.node_type, 'UnaryOpNode')
        ols_error('expecting Variable or UnaryOp on LHS', line);
    else
        if strcmp(ast{i}.AST.arg1.node_type, 'VariableNode') ...
                && (ast{i}.AST.arg1.lag ~= 0 ...
                || ~any(strcmp(ds.name, ast{i}.AST.arg1.name)))
            ols_error('the lhs of the equation must be an Variable or UnaryOp with lag == 0 that exists in the dataset', line);
        end
        if strcmp(ast{i}.AST.arg1.node_type, 'UnaryOpNode') ...
                && (ast{i}.AST.arg1.arg.lag ~= 0 ...
                || ~any(strcmp(ds.name, ast{i}.AST.arg1.arg.name)))
            ols_error('the lhs of the equation must be an Variable or UnaryOp with lag == 0 that exists in the dataset', line);
        end
    end

    % Set LHS (Y)
    lhssub = dseries();
    Y = evalNode(ds, ast{i}.AST.arg1, jsonmodel{i}.line, dseries());

    % Set RHS (X)
    plus_node = ast{i}.AST.arg2;
    last_node_to_parse = [];
    residual = '';
    X = dseries();
    while ~isempty(plus_node) || ~isempty(last_node_to_parse)
        Xtmp = dseries();
        if isempty(last_node_to_parse)
            [plus_node, node_to_parse, last_node_to_parse] = findNextplus_node(plus_node, jsonmodel{i}.line);
        else
            node_to_parse = last_node_to_parse;
            last_node_to_parse = [];
        end
        if strcmp(node_to_parse.node_type, 'VariableNode') || strcmp(node_to_parse.node_type, 'UnaryOpNode')
            if strcmp(node_to_parse.node_type, 'VariableNode')
                if strcmp(node_to_parse.type, 'parameter')
                    % Intercept
                    Xtmp = dseries(ones(ds.nobs, 1), ds.firstdate, node_to_parse.name);
                elseif strcmp(node_to_parse.type, 'exogenous') && ~any(strcmp(ds.name, node_to_parse.name))
                    % Residual if not contained in ds
                    if isempty(residual)
                        residual = node_to_parse.name;
                    else
                        ols_error(['only one residual allowed per equation; encountered ' residual ' & ' node_to_parse.name], jsonmodel{i}.line);
                    end
                elseif strcmp(node_to_parse.type, 'endogenous') ...
                        || (strcmp(node_to_parse.type, 'exogenous') && any(strcmp(ds.name, node_to_parse.name)))
                    % Subtract VariableNode from LHS
                    % NB: treat exogenous that exist in ds as endogenous
                    lhssub = lhssub + evalNode(ds, node_to_parse, jsonmodel{i}.line, dseries());
                else
                    ols_error('unexpected variable type found', jsonmodel{i}.line);
                end
            else
                % Subtract UnaryOpNode from LHS
                % NB: treat exogenous that exist in ds as endogenous
                lhssub = lhssub + evalNode(ds, node_to_parse, jsonmodel{i}.line, dseries());
            end
        elseif strcmp(node_to_parse.node_type, 'BinaryOpNode') && strcmp(node_to_parse.op, '*')
            % Parse param_expr * endog_expr
            Xtmp = parseTimesNode(ds, node_to_parse, jsonmodel{i}.line);
            if length(Xtmp.name) > 1
                % Handle constraits
                % Look through Xtmp names for constant
                % if found, subtract from LHS
                to_remove = [];
                for j = 1:length(Xtmp.name)
                    if ~isnan(str2double(Xtmp.name{j}))
                        lhssub = lhssub + str2double(Xtmp.name{j}) * Xtmp{j};
                        to_remove = [to_remove j];
                    end
                    % Multiply by -1 now so that it can be added together below
                    % Otherwise, it would matter which was encountered first,
                    % a parameter on its own or a linear constraint
                    Xtmp.(Xtmp.name{j}) = -1 * Xtmp{j};
                end
                for j = length(to_remove):-1:1
                    Xtmp = Xtmp.remove(Xtmp.name{j});
                end
            end
        else
            ols_error('didn''t expect to arrive here', jsonmodel{i}.line);
        end
        if ~isempty(Xtmp)
            to_remove = [];
            for j = 1:length(Xtmp.name)
                % Handle constraits
                idx = find(strcmp(X.name, Xtmp.name{j}));
                if ~isempty(idx)
                    X.(X.name{idx}) = X{idx} + Xtmp{j};
                    to_remove = [to_remove j];
                end
            end
            for j = length(to_remove):-1:1
                Xtmp = Xtmp.remove(Xtmp.name{j});
            end
            X = [X Xtmp];
        end
    end

    clear residual Xtmp plus_node last_node_to_parse to_remove

    Y = Y - lhssub;

    %% Check to see if called from Gibbs
    st = dbstack(1);
    varargout = cell(1, 1);
    called_from_olsgibbs = false;
    if strcmp(st(1).name, 'olsgibbs')
        varargout = cell(1, 4);
        called_from_olsgibbs = true;
    end

    %% Set dates
    fp = max(Y.firstobservedperiod, X.firstobservedperiod);
    lp = min(Y.lastobservedperiod, X.lastobservedperiod);
    if ~isempty(lhssub)
        fp = max(fp, lhssub.firstobservedperiod);
        lp = min(lp, lhssub.lastobservedperiod);
    end
    if isfield(jsonmodel{i}, 'tags') ...
            && isfield(jsonmodel{i}.tags, 'sample') ...
            && ~isempty(jsonmodel{i}.tags.sample)
        colon_idx = strfind(jsonmodel{i}.tags.sample, ':');
        fsd = dates(jsonmodel{i}.tags.sample(1:colon_idx-1));
        lsd = dates(jsonmodel{i}.tags.sample(colon_idx+1:end));
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
    if ~isempty(lhssub)
        lhssub = lhssub(fp:lp);
    end
    pnames = X.name;
    X = X(fp:lp).data;

    [nobs, nvars] = size(X);
    if called_from_olsgibbs
        varargout{1} = nobs;
        varargout{2} = pnames;
        varargout{3} = X;
        varargout{4} = Y.data;
        varargout{5} = jsonmodel;
        varargout{6} = fp;
        varargout{7} = lp;
        return
    end

    if isfield(jsonmodel{i}, 'tags') && ...
            isfield(jsonmodel{i}.tags, 'name')
        tag = jsonmodel{i}.tags.('name');
    else
        tag = ['eq_line_no_' num2str(jsonmodel{i}.line)];
    end

    %% Estimation
    % From LeSage, James P. "Applied Econometrics using MATLAB"
    oo_.ols.(tag).dof = nobs - nvars;

    % Estimated Parameters
    [q, r] = qr(X, 0);
    xpxi = (r'*r)\eye(nvars);
    oo_.ols.(tag).beta = r\(q'*Y.data);
    oo_.ols.(tag).param_idxs = zeros(length(pnames), 1);
    for j = 1:length(pnames)
        if ~strcmp(pnames{j}, 'intercept')
            oo_.ols.(tag).param_idxs(j) = find(strcmp(M_.param_names, pnames{j}));
            M_.params(oo_.ols.(tag).param_idxs(j)) = oo_.ols.(tag).beta(j);
        end
    end

    % Yhat
    idx = 0;
    yhatname = [tag '_FIT'];
    if ~isempty(fitted_names_dict)
        idx = strcmp(fitted_names_dict(:,1), tag);
        if any(idx)
            yhatname = fitted_names_dict{idx, 2};
        end
    end
    oo_.ols.(tag).Yhat = dseries(X*oo_.ols.(tag).beta, fp, yhatname);

    % Residuals
    oo_.ols.(tag).resid = Y - oo_.ols.(tag).Yhat;

    % Correct Yhat reported back to user
    Y = Y + lhssub;
    oo_.ols.(tag).Yhat = oo_.ols.(tag).Yhat + lhssub;

    % Apply correcting function for Yhat if it was passed
    if any(idx) ...
            && length(fitted_names_dict(idx, :)) == 3 ...
            && ~isempty(fitted_names_dict{idx, 3})
        oo_.ols.(tag).Yhat = ...
            feval(fitted_names_dict{idx, 3}, oo_.ols.(tag).Yhat);
    end
    ds = [ds oo_.ols.(tag).Yhat];

    %% Calculate statistics
    % Estimate for sigma^2
    SS_res = oo_.ols.(tag).resid.data'*oo_.ols.(tag).resid.data;
    oo_.ols.(tag).s2 = SS_res/oo_.ols.(tag).dof;

    % R^2
    ym = Y.data - mean(Y);
    SS_tot = ym'*ym;
    oo_.ols.(tag).R2 = 1 - SS_res/SS_tot;

    % Adjusted R^2
    oo_.ols.(tag).adjR2 = oo_.ols.(tag).R2 - (1 - oo_.ols.(tag).R2)*(nvars-1)/(oo_.ols.(tag).dof);

    % Durbin-Watson
    ediff = oo_.ols.(tag).resid.data(2:nobs) - oo_.ols.(tag).resid.data(1:nobs-1);
    oo_.ols.(tag).dw = (ediff'*ediff)/SS_res;

    % Standard Error
    oo_.ols.(tag).stderr = sqrt(oo_.ols.(tag).s2*diag(xpxi));

    % T-Stat
    oo_.ols.(tag).tstat = oo_.ols.(tag).beta./oo_.ols.(tag).stderr;

    %% Print Output
    if ~options_.noprint
        if nargin == 3
            title = ['OLS Estimation of equation ''' tag ''' [name = ''' tag ''']'];
        else
            title = ['OLS Estimation of equation ''' tag ''''];
        end

        preamble = {sprintf('Dependent Variable: %s', jsonmodel{i}.lhs), ...
            sprintf('No. Independent Variables: %d', nvars), ...
            sprintf('Observations: %d from %s to %s\n', nobs, fp.char, lp.char)};

        afterward = {sprintf('R^2: %f', oo_.ols.(tag).R2), ...
            sprintf('R^2 Adjusted: %f', oo_.ols.(tag).adjR2), ...
            sprintf('s^2: %f', oo_.ols.(tag).s2), ...
            sprintf('Durbin-Watson: %f', oo_.ols.(tag).dw)};

        dyn_table(title, preamble, afterward, pnames, ...
            {'Estimates','t-statistic','Std. Error'}, 4, ...
            [oo_.ols.(tag).beta oo_.ols.(tag).tstat oo_.ols.(tag).stderr]);
    end

    if ~called_from_olsgibbs
        varargout{1} = ds;
    end
end
end

%% Helper Functions
function ols_error(msg, line)
error(['ERROR encountered in `dyn_ols` line ' num2str(line) ': ' msg])
end

function [next_plus_node, node_to_parse, last_node_to_parse] = findNextplus_node(plus_node, line)
% Given an additive entry in the AST, find the next additive entry
% (next_plus_node). Also find the node that will be parsed into
% parameter*endogenous||param||exog|endog (node_to_parse).
% Function used for moving through the AST.
if ~(strcmp(plus_node.node_type, 'BinaryOpNode') && strcmp(plus_node.op, '+'))
    ols_error('pairs of nodes must be separated additively', line);
end
next_plus_node = [];
last_node_to_parse = [];
if strcmp(plus_node.arg1.node_type, 'BinaryOpNode') && strcmp(plus_node.arg1.op, '+')
    next_plus_node = plus_node.arg1;
    node_to_parse = getOlsNode(plus_node.arg2, line);
elseif strcmp(plus_node.arg2.node_type, 'BinaryOpNode') && strcmp(plus_node.arg2.op, '+')
    next_plus_node = plus_node.arg2;
    node_to_parse = getOlsNode(plus_node.arg1, line);
else
    node_to_parse = getOlsNode(plus_node.arg1, line);
    last_node_to_parse = getOlsNode(plus_node.arg2, line);
end
end

function node_to_parse = getOlsNode(node, line)
if ~(strcmp(node.node_type, 'BinaryOpNode') && strcmp(node.op, '*')) ...
        && ~strcmp(node.node_type, 'VariableNode') ...
        && ~strcmp(node.node_type, 'UnaryOpNode')
    ols_error('couldn''t find node to parse', line);
end
node_to_parse = node;
end

function X = parseTimesNode(ds, node, line)
% Separate the parameter expression from the endogenous expression
assert(strcmp(node.node_type, 'BinaryOpNode') && strcmp(node.op, '*'))
[param, X] = parseTimesNodeHelper(ds, node.arg1, line, {}, dseries());
[param, X] = parseTimesNodeHelper(ds, node.arg2, line, param, X);
X = X.rename(param{1});
for ii = 2:length(param)
    X = [X dseries(X{1}.data, X{1}.firstdate, param{ii})];
end
end

function [param, X] = parseTimesNodeHelper(ds, node, line, param, X)
if isOlsParamExpr(node, line)
    param = assignParam(param, node, line);
elseif isOlsVarExpr(node, line)
    if isempty(X)
        X = evalNode(ds, node, line, X);
    else
        ols_error(['got endog * endog' node.name ' (' node.type ')'], line);
    end
else
    ols_error('unexpected expression', line);
end
end

function param = assignParam(param, node, line)
if ~isempty(param)
    ols_error(['got param * param' node.name ' (' node.type ')'], line);
end
param = assignParamHelper(param, node, line);
end

function param = assignParamHelper(param, node, line)
if strcmp(node.node_type, 'NumConstNode')
    param{end+1} = num2str(node.value);
elseif strcmp(node.node_type, 'VariableNode')
    param{end+1} = node.name;
elseif strcmp(node.node_type, 'BinaryOpNode')
    if ~strcmp(node.op, '-')
        ols_error(['got unexpected parameter op ' node.op], line);
    end
    param = assignParamHelper(param, node.arg1, line);
    param = assignParamHelper(param, node.arg2, line);
else
    ols_error(['got unexpected node (' node.type ')'], line);
end
end

function tf = isOlsVar(node)
if strcmp(node.node_type, 'VariableNode') ...
        && (strcmp(node.type, 'endogenous') ...
        || (strcmp(node.type, 'exogenous') && any(strcmp(ds.name, node.name))))
    tf = true;
elseif strcmp(node.node_type, 'UnaryOpNode') ...
        && (strcmp(node.arg.type, 'endogenous') ...
        || (strcmp(node.arg.type, 'exogenous') && any(strcmp(ds.name, node.arg.name))))
    tf = true;
else
    tf = false;
end
end

function tf = isOlsVarExpr(node, line)
if strcmp(node.node_type, 'VariableNode') || strcmp(node.node_type, 'UnaryOpNode')
    tf = isOlsVar(node);
elseif strcmp(node.node_type, 'BinaryOpNode')
    tf = isOlsVarExpr(node.arg1, line) || isOlsVarExpr(node.arg2, line);
else
    ols_error(['got unexpected type ' node.node_type], line);
end
end

function X = evalNode(ds, node, line, X)
if strcmp(node.node_type, 'NumConstNode')
    X = dseries(str2double(node.value), ds.dates, 'const');
elseif strcmp(node.node_type, 'VariableNode')
    if ~(strcmp(node.type, 'endogenous') ...
            || (strcmp(node.type, 'exogenous') && any(strcmp(ds.name, node.name))))
        ols_error(['got unexpected type ' node.name ': ' node.type], line);
    end
    X = ds.(node.name)(node.lag);
elseif strcmp(node.node_type, 'UnaryOpNode')
    Xtmp = evalNode(ds, node.arg, line, X);
    % Only works if dseries supports . notation for unary op (true for log/diff)
    % Otherwise, use: X = eval([node.op '(Xtmp)']);
    X = Xtmp.(node.op);
elseif strcmp(node.node_type, 'BinaryOpNode')
    Xtmp1 = evalNode(ds, node.arg1, line, X);
    Xtmp2 = evalNode(ds, node.arg2, line, X);
    X = X + eval(['Xtmp1 ' node.op ' Xtmp2']);
else
    ols_error(['got unexpected node type ' node.node_type], line);
end
end

function tf = isOlsParam(node)
if strcmp(node.node_type, 'VariableNode') && strcmp(node.type, 'parameter')
    tf = true;
else
    tf = false;
end
end

function tf = isOlsParamExpr(node, line)
if strcmp(node.node_type, 'NumConstNode')
    tf = true;
elseif strcmp(node.node_type, 'VariableNode')
    tf = isOlsParam(node);
elseif strcmp(node.node_type, 'UnaryOpNode')
    tf = false;
elseif strcmp(node.node_type, 'BinaryOpNode')
    tf = isOlsParamExpr(node.arg1) && isOlsParamExpr(node.arg2);
    if tf && ~strcmp(node.op, '-')
        ols_error(['got unexpected op ' node.op], line);
    end
else
    ols_error(['got unexpected type ' node.node_type], line);
end
end
