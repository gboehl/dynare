function aggregate(ofile, dynopt, rootfolder, varargin)

% Agregates cherry-picked models.

% Copyright (C) 2019-2021 Dynare Team
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

MAX_NUMBER_OF_ELEMENTS = 10000;

warning off MATLAB:subscripting:noSubscriptsSpecified

if ~isempty(dynopt)
    % Should be a list of options for the preprocessor in a cell
    % array.
    firstline = '// --+ options:';
    if iscell(dynopt)
        for i = 1:length(dynopt)
            firstline = sprintf('%s %s', firstline, dynopt{i});
        end
        firstline = sprintf('%s %s', firstline, '+--');
    else
        error('Second argument has to be a cell array (list of options for dynare preprocessor).')
    end
else
    firstline = '';
end


% Get parameters.
for i=1:length(varargin)
    fid = fopen(sprintf('%s/parameters.inc', varargin{i}));
    if fid<0
        % No parameters in the cherrypicked (sub)model, go to the
        % next cherrypicked model.
        continue
    end
    statement = fgetl(fid);
    if exist('plist', 'var')
        plist = union(plist, strsplit(statement, {'parameters', ' ', ';'}));
    else
        plist = strsplit(statement, {'parameters', ' ', ';'});
    end
    plist(cellfun(@(x) all(isempty(x)), plist)) = [];
    fclose(fid);
end

% Get equations
eqlist = cell(MAX_NUMBER_OF_ELEMENTS, 4);
tagnum = 1;
eqnum = 0;
for i=1:length(varargin)
    % Store all non-empty lines of model.inc in the “model” cell-array
    fid = fopen(sprintf('%s/model.inc', varargin{i}));
    model = {};
    while ~feof(fid)
        line = fgetl(fid);
        if ~isempty(line)
            model{end+1} = line;
        end
    end
    fclose(fid);

    eqtag = false;
    for j=1:length(model)
        if isequationtag(model{j})
            if eqtag
                error('An equation tag must be followed by an equation.')
            end
            % Ensure that the equation tag name matches the LHS variable.
            eqtagname = regexp(model{j}, 'name=''(\w*)''', 'match');
            [lhs, ~] = getequation(model{j+1});
            endovar = getendovar(lhs);
            eqtagname_ = strcat('name=''', endovar{1}, '''');
            if ~isempty(eqtagname)
                if ~isequal(eqtagname{1}, eqtagname_)
                    model{j} = strrep(model{j}, eqtagname{1}, eqtagname_);
                end
            else
                model{j} = eqtagname_;
            end
            % Add equation tag with block name.
            if ~isempty(rootfolder)
                model{j} = strcat('[blockname=''',  getblockname(varargin{i}, rootfolder), ''',', model{j}(2:end));
            end
            eqlist{tagnum,4} = model{j};
            eqtag = true;
        else
            eqnum = eqnum+1;
            [lhs, rhs] = getequation(model{j});
            endovar = getendovar(lhs);
            eqlist{eqnum,1} = endovar{1};
            eqlist{eqnum,2} = lhs;
            eqlist{eqnum,3} = rhs;
            eqtag = false;
            tagnum = tagnum+1;
        end
    end
end
eqlist = eqlist(1:eqnum,:);
if isoctave && octave_ver_less_than('6')
    [~, idx] = unique_stable(eqlist(:,1));
else
    [~, idx] = unique(eqlist(:,1), 'stable');
end
eqlist = eqlist(idx, :);

% Get endogenous variables.
elist = cell(MAX_NUMBER_OF_ELEMENTS, 2);
enum = 0;
for i=1:length(varargin)
    fid = fopen(sprintf('%s/endogenous.inc', varargin{i}));
    cline = fgetl(fid);
    while ischar(cline)
        if ~isequal(cline, 'var')
            enum = enum+1;
            cline = regexprep(cline, '\t', '');
            cline = regexprep(cline, ';', '');
            [v, t] = getvarandtag(cline);
            elist(enum,1) = {v};
            elist(enum,2) = {t};
        end
        cline = fgetl(fid);
    end
    fclose(fid);
end
elist = elist(1:enum,:);
if isoctave && octave_ver_less_than('6')
    [~, idx] = unique_stable(elist(:,1));
else
    [~, idx] = unique(elist(:,1), 'stable');
end
elist = elist(idx,:);

% Get exogenous variables.
xlist = cell(MAX_NUMBER_OF_ELEMENTS, 2);
xnum = 0;
for i=1:length(varargin)
    fid = fopen(sprintf('%s/exogenous.inc', varargin{i}));
    if fid<0
        % No exogenous variables in the cherrypicked (sub)model, go to the
        % next cherrypicked model.
        continue
    end
    cline = fgetl(fid);
    while ischar(cline)
        if ~isequal(cline, 'varexo')
            xnum = xnum+1;
            cline = regexprep(cline, '\t', '');
            cline = regexprep(cline, ';', '');
            [v, t] = getvarandtag(cline);
            xlist(xnum,1) = {v};
            xlist(xnum,2) = {t};
        end
        cline = fgetl(fid);
    end
    fclose(fid);
end
xlist = xlist(1:xnum,:);
if isoctave && octave_ver_less_than('6')
    [~, idx] = unique_stable(xlist(:,1));
else
    [~, idx] = unique(xlist(:,1), 'stable');
end
xlist = xlist(idx,:);

% Get parameter values.
calibration = '';
for i=1:length(varargin)
    fid = fopen(sprintf('%s/parameter-values.inc', varargin{i}));
    if fid<0
        % No calibrations in the cherrypicked (sub)model, go to the
        % next cherrypicked model.
        continue
    end
    cline = fgetl(fid);
    while ischar(cline)
        calibration = sprintf('%s\n%s', calibration, cline);
        cline = fgetl(fid);
    end
    fclose(fid);
end

% Move the endogenous variables which are not LHS of an equation
% into the set of exogenous variables.
[~, i1] = intersect(elist(:,1), eqlist(:,1));
if ~isequal(length(i1),rows(eqlist))
    error('Something is wrong with the endogenous variables.')
end
i2 = setdiff(1:rows(elist), i1);
xlist = [xlist; elist(i2,:)];
[~,idx] = unique(xlist(:,1));          % Ensure that the exogenous variable names are unique.
xlist = [xlist(idx,1) xlist(idx,2)];   % We do not test that the tags are the same.
elist = elist(i1,:);

% Remove endogenous variables from list of exogenous variables (if any).
xlist1 = xlist(:,1);
xlist2 = xlist(:,2);
[xlist1, id] = setdiff(xlist1, elist(:,1));
xlist2 = xlist2(id);
xlist = [xlist1, xlist2];

% Print all cherry-picked models in one mod-file.
[filepath, filename, fileext] = fileparts(ofile);
if ~isempty(filepath) && ~exist(filepath, 'dir')
    mkdir(filepath);
end
if isempty(filepath)
    fid = fopen(sprintf('%s%s', filename, fileext), 'w');
else
    fid = fopen(sprintf('%s%s%s%s', filepath, filesep(), filename, fileext), 'w');
end
if ~isempty(firstline)
    fprintf(fid, '%s\n\n', firstline);
end
fprintf(fid, 'var\n');
for i=1:rows(elist)
    if size(elist,2)==1 || isempty(elist{i,2})
        fprintf(fid, '\t%s\n', elist{i,1});
    else
        fprintf(fid, '\t%s %s\n', elist{i,1}, elist{i,2});
    end
end
if ~isempty(plist)
    fprintf(fid, ';\n\n');
    fprintf(fid, 'parameters\n');
    for i=1:length(plist)
        fprintf(fid, '\t%s\n', plist{i});
    end
    fprintf(fid, ';\n\n');
    fprintf(fid, calibration);
end
if ~isempty(xlist)
    fprintf(fid, '\n\n');
    fprintf(fid, 'varexo\n');
    for i=1:rows(xlist)
        if size(xlist,2)==1 || isempty(xlist{i,2})
            fprintf(fid, '\t%s\n', xlist{i,1});
        else
            fprintf(fid, '\t%s %s\n', xlist{i,1}, xlist{i,2});
        end
    end
end
fprintf(fid, ';\n\n');
fprintf(fid, 'model;\n\n');
for i=1:rows(eqlist)
    if isempty(eqlist{i,4})
        fprintf(fid, '\t%s = %s;\n\n', eqlist{i,2}, eqlist{i,3});
    else
        fprintf(fid, '\t%s\n', eqlist{i,4});
        fprintf(fid, '\t%s = %s;\n\n', eqlist{i,2}, eqlist{i,3});
    end
end
fprintf(fid, 'end;');
fclose(fid);

warning on MATLAB:subscripting:noSubscriptsSpecified


function b = isequationtag(str)
    b = true;
    if isempty(regexp(str, '\[.*\]','once'))
        b = false;
    end

function [lhs, rhs] = getequation(str)
    terms = strsplit(str, {'=',';'});
    terms(cellfun(@(x) all(isempty(x)), terms)) = [];
    terms(1) = {strrep(terms{1}, ' ', '')};
    lhs = regexp(terms{1}, '^(diff\([\-]*(log|diff)\([\-\+\*\/\w]*\)\)|(log|diff)\([\(\-\+\*\/\)\w]*\)|\w*)', 'match');
    if ~isempty(lhs)
        lhs = lhs{1};
        if isequal(lhs, 'log')
           error('Malformed equation: log of log or diff are not allowed.')
        end
        rhs = terms{2};
    else
        error('Malformed equation.')
    end

function v = getendovar(lhs)
    v = strsplit(lhs, {'diff','log','(',')', '+', '-', '*', '/'});
    v(cellfun(@(x) all(isempty(x)), v)) = [];
    if length(v)>1
        error('Malformed equation: no more than one endogenous variable can be used on the LHS.')
    end

function [v, t] = getvarandtag(str)
    tmp = regexp(str, '(?<name>\w+)\s*(?<tag>\(.*\))', 'names');
    if isempty(tmp)
        tmp = regexp(str, '(?<name>\w+)\s*', 'names');
        v = tmp.name;
        t = '';
    else
        v = tmp.name;
        t = tmp.tag;
    end

function blkname = getblockname(str, ROOT_FOLDER)
    str = strrep(str, '/', filesep());
    str = strrep(str, [ROOT_FOLDER filesep() 'blocks' filesep()], '');
    idx = strfind(str, filesep());
    blkname = str(1:idx(1)-1);
