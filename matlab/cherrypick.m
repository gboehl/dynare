function json = cherrypick(infile, outfold, eqtags, noresids, json)

% Extract some equations in infile (mod file used for estimation)
% and write them in outfile (mod file used for simulation).
%
% INPUTS
% - infile        [string]    Name of the mod file where all the equations used for estimation are available.
% - outfold       [string]    Name of the folder where the generated files are saveda subset of the equations is to be printed.
% - eqtags        [cell]      Equation tags of the selected equations.
% - noresids      [logical]   Removes estimation residuals (not to be used in simulation) if true.
% - json          [char]      Content of a JSON file.
%
% OUTPUTS
% - json          [char]      Content of a JSON file.
%
% SPECIAL REQUIREMENTS
% It is expected that the file infile.mod has already been run, and
% that the associated JSON output is available.

% Copyright © 2019-2021 Dynare Team
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

global M_

% Set default value
if nargin<4 || isempty(noresids)
    noresids = true;
end

% Delete outfold subdirectory if it already exists
if exist(outfold, 'dir')
    rmdir(outfold, 's');
end

% Create the subdirectoty where the generated files will be saved.
mkdir(outfold);

% Check that infile.mod and the related JSON output exist.
if ~exist(sprintf('%s.mod', infile), 'file')
    error('Cannot find %s.mod.', infile)
end
if ~exist(sprintf('%s/model/json', infile), 'dir')
    error('Cannot find %s/model/json folder. Did you run %s.mod with the json option?', infile, infile);
end

% Check if some variables have to be renamed.
rename = M_.equations_tags(strcmp('rename',M_.equations_tags(:,2)),[1,3]);
isrename = ~isempty(rename);

if nargin<5
    % Load json file (original mod file)
    json = loadjson_(sprintf('%s/model/json/modfile-original.json', M_.dname));
end

% Create a new file.
fid = fopen(sprintf('%s/model.inc', outfold), 'w');

plist = {};
elist = {};
xlist = {};

for i=1:length(eqtags)
    rhs = [];
    lhs = [];
    % Get equation number.
    eqnum = get_equation_number_by_tag(eqtags{i}, M_);
    % Get the original equation.
    [LHS, RHS] = get_lhs_and_rhs(eqtags{i}, M_, true, json);
    % Get the parameters, endogenous and exogenous variables in the current equation.
    [pnames, ~, xnames] = get_variables_and_parameters_in_equation(LHS, RHS, M_);
    lhs_expression = LHS;
    LHS = get_variables_and_parameters_in_expression(LHS);
    enames = LHS;
    if length(LHS)>1
        error('Expressions with more than one variable on the LHS are not allowed.')
    end
    LHS = LHS{1};
    if isrename
        [variable_has_to_be_renamed, id] = ismember(eqnum, [rename{:,1}]);
        if variable_has_to_be_renamed
            TMP = strsplit(rename{id,2}, ',');
            for j=1:length(TMP)
                tmp = strsplit(TMP{j}, '->');
                lhs_expression = exactstrrep(lhs_expression, tmp{1}, tmp{2});
                RHS = exactstrrep(RHS, tmp{1}, tmp{2});
                rep = strcmp(tmp{1}, enames);
                if any(rep)
                    enames(rep) = tmp(2);
                end
                rep = strcmp(tmp{1}, xnames);
                if any(rep)
                    xnames(rep) = tmp(2);
                end
            end
        end
    end
    % Remove residual from equation if required.
    if noresids
        exogenous_variables_to_be_removed = ~ismember(xnames, M_.simulation_exo_names);
        if any(exogenous_variables_to_be_removed)
            switch sum(exogenous_variables_to_be_removed)
              case 1
                RHS = regexprep(RHS, sprintf('\\ *\\+\\ *%s', xnames{exogenous_variables_to_be_removed}), '');
                RHS = regexprep(RHS, sprintf('%s', xnames{exogenous_variables_to_be_removed}), '');
              case 0
                % Nothing to do.
              otherwise
                error('Cannot remove more than one exogenous variable in an equation (%s).', eqtags{i})
            end
            xnames = setdiff(xnames, xnames{exogenous_variables_to_be_removed});
        end
    end
    % Unroll expectation terms if any.
    isvar = regexp(RHS, 'var_expectation\(model_name = (?<name>\w+)\)', 'names');
    ispac = regexp(RHS, 'pac_expectation\(model_name = (?<name>\w+)\)', 'names');
    if ~isempty(isvar)
        rhs = write_expectations(eqtags{i}, isvar.name, 'var');
        lhs = sprintf('%s_VE', eqtags{i});
        RHS = strrep(RHS, sprintf('var_expectation(model_name = %s)', isvar.name), lhs);
    end
    if ~isempty(ispac)
        [rhs, growthneutralitycorrection] = write_expectations(eqtags{i}, ispac.name, 'pac');
        if ~isempty(rhs)
            lhs = sprintf('%s_PE', eqtags{i});
            if isempty(growthneutralitycorrection)
                RHS = strrep(RHS, sprintf('pac_expectation(model_name = %s)', ispac.name), lhs);
            else
                RHS = strrep(RHS, sprintf('pac_expectation(model_name = %s)', ispac.name), sprintf('%s+%s', lhs, growthneutralitycorrection));
            end
        else
            % MCE version of the PAC equation.
            [rhs, growthneutralitycorrection] = write_pac_mce_expectations(eqtags{i}, ispac.name);
            lhs = sprintf('%s_Z', eqtags{i});
            if isempty(growthneutralitycorrection)
                RHS = strrep(RHS, sprintf('pac_expectation(model_name = %s)', ispac.name), lhs);
            else
                RHS = strrep(RHS, sprintf('pac_expectation(model_name = %s)', ispac.name), sprintf('%s+%s', lhs, growthneutralitycorrection));
            end
        end
    end
    % Print equation for unrolled PAC/VAR-expectation and update
    % list of parameters and endogenous variables (if any).
    if ~isempty(rhs)
        % Note that the call to get_variables_and_parameters_in_equation()
        % will not return the lhs variable in expectation_enames since
        % the name is created on the fly and is not a  member of M_.endo_names.
        expectation_pnames = get_variables_and_parameters_in_equation('', rhs, M_);
        expectation_enames = get_variables_and_parameters_in_expression(lhs);
        expectation_xnames = get_variables_and_parameters_in_expression(rhs);
        pnames = union(pnames, expectation_pnames);
        xnames = union(xnames, setdiff(expectation_xnames, expectation_pnames));
        enames = union(enames, expectation_enames);
        fprintf(fid, '[name=''%s'']\n', lhs);
        fprintf(fid, '%s = %s;\n\n', lhs, rhs);
    else
        pRHS = get_variables_and_parameters_in_equation('', RHS, M_);
        xRHS = get_variables_and_parameters_in_expression(RHS);
        xnames = union(xnames, setdiff(xRHS, pRHS));
        pnames = union(pnames, pRHS);
    end
    % Update pnames, enames and xnames if PAC with growth neutrality correction.
    if ~isempty(ispac) && ~isempty(growthneutralitycorrection)
        [growthneutralitycorrection_pnames, ...
             growthneutralitycorrection_enames, ...
             growthneutralitycorrection_xnames] = get_variables_and_parameters_in_equation('', growthneutralitycorrection, M_);
        if ~isempty(growthneutralitycorrection_pnames)
            pnames = union(pnames, growthneutralitycorrection_pnames);
        end
        if ~isempty(growthneutralitycorrection_enames)
            xnames = union(xnames, growthneutralitycorrection_enames);
        end
        if ~isempty(growthneutralitycorrection_xnames)
            xnames = union(xnames, growthneutralitycorrection_xnames);
        end
    end
    % Print tags
    if iscell(json.model)
        tfields = fieldnames(json.model{eqnum}.tags);
        tags = sprintf('%s=''%s''', tfields{1}, json.model{eqnum}.tags.(tfields{1}));
        for j=2:length(tfields)
            if ~isempty(json.model{eqnum}.tags.(tfields{j}))
                tags = sprintf('%s, %s=''%s''', tags, tfields{j}, json.model{eqnum}.tags.(tfields{j}));
            end
        end
    else
        tfields = fieldnames(json.model.tags);
        tags = sprintf('%s=''%s''', tfields{1}, json.model.tags.(tfields{1}));
        for j=2:length(tfields)
            if ~isempty(json.model.tags.(tfields{j}))
                tags = sprintf('%s, %s=''%s''', tags, tfields{j}, json.model.tags.(tfields{j}));
            end
        end
    end
    fprintf(fid, '[%s]\n', tags);
    % Print equation.
    fprintf(fid, '%s = %s;\n\n', lhs_expression, RHS);
    % Update lists of parameters, endogenous variables and exogenous variables.
    plist = union(plist, pnames);
    elist = union(elist, enames);
    xlist = union(xlist, xnames);
end
fclose(fid);

% Export parameters
if ~isempty(plist)
    fid = fopen(sprintf('%s/parameters.inc', outfold), 'w');
    fprintf(fid, 'parameters %s;', sprintf('%s ', plist{:}));
    fclose(fid);
end

% Export endogegnous variables
fid = fopen(sprintf('%s/endogenous.inc', outfold), 'w');
printlistofvariables(fid, 'endo', elist, M_, elist);
fclose(fid);

% Export exogenous variables
if ~isempty(xlist)
    fid = fopen(sprintf('%s/exogenous.inc', outfold), 'w');
    printlistofvariables(fid, 'exo', xlist, M_, xlist);
    fclose(fid);
end

% Export parameter values
if ~isempty(plist)
    fid = fopen(sprintf('%s/parameter-values.inc', outfold), 'w');
    for i=1:length(plist)
        id = strcmp(plist{i}, M_.param_names);
        if any(id)
            fprintf(fid, '%s = %s;\n', plist{i}, num2str(M_.params(id), 16));
        end
    end
    fclose(fid);
end

function printlistofvariables(fid, kind, list, DynareModel, vappend)
    if isfield(DynareModel, sprintf('%s_partitions', kind))
        % Some endogenous variables are tagged.
        switch kind
          case 'exo'
            tfields = fieldnames(DynareModel.exo_partitions);
            vlist = 'varexo';
            vnames = DynareModel.exo_names;
            partitions = DynareModel.exo_partitions;
          case 'endo'
            tfields = fieldnames(DynareModel.endo_partitions);
            vlist = 'var';
            vnames = DynareModel.endo_names(1:DynareModel.orig_endo_nbr);
            partitions = DynareModel.endo_partitions;
          otherwise
            error('Illegal value for second input argument.')
        end
        for i = 1:length(list)
            id = strmatch(list{i}, vnames, 'exact');
            if ~isempty(id)
                tags = '';
                for j=1:length(tfields)
                    if ~isempty(partitions.(tfields{j}){id})
                        tags = sprintf('%s, %s=''%s''', tags, tfields{j}, partitions.(tfields{j}){id});
                    end
                end
                if ~isempty(tags)
                    tags = sprintf('(%s)', tags(3:end));
                end
            elseif ~isempty(strmatch(list{i}, vappend, 'exact'))
                % Nothing to do, this variable was renamed by cherrypick
                tags = '';
            else
                if isequal(kind, 'endo') &&  (isequal(list{i}(end-2:end), '_PE') || isequal(list{i}(end-2:end), '_VE'))
                    if isequal(list{i}(end-2:end), '_PE')
                        tags = sprintf('(expectation_kind=''%s'')', 'pac');
                    else
                        tags = sprintf('(expectation_kind=''%s'')', 'var');
                    end
                else
                    error('Unknown variable.')
                end
            end
            if isempty(tags)
                vlist = sprintf('%s\n\t%s', vlist, list{i});
            else
                vlist = sprintf('%s\n\t%s %s', vlist, list{i}, tags);
            end
        end
        fprintf(fid, '%s;', vlist);
    else
        switch kind
          case 'exo'
            vlist = 'varexo';
          case 'endo'
            vlist = 'var';
          otherwise
            error('Illegal value for second input argument.')
        end
        for i=1:length(list)
            vlist = sprintf('%s\n\t%s', vlist, list{i});
        end
        fprintf(fid, '%s;', vlist);
    end
